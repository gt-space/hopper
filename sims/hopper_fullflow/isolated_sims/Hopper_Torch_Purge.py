"""Run just the torch cold-gas purge, without a propulsion-system solve.

python isolated_sims/Hopper_Torch_Purge.py [--config case.json] [--out folder]
Edit the entire timeline in torch.standalone_purge.sequence in hopper_config.json.
"""
import argparse
import copy
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys
import time

if not __package__:
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))

from Backend.Hopper_Config import load_config
from Backend.Hopper_Torch_Purge import TorchPurgeModel


def write_csv(path, rows):
    with path.open('w',newline='',encoding='utf-8') as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def archive_config(config, source_path, output):
    """Keep source text and a replayable config with referenced input assets."""
    (output/'input_config.json').write_text(source_path.read_text(encoding='utf-8-sig'),encoding='utf-8')
    case=copy.deepcopy(config)
    assets=output/'inputs'
    assets.mkdir()
    def archive(relative, name):
        source=(source_path.parent/relative).resolve()
        destination=assets/name
        shutil.copy2(source,destination)
        return destination.relative_to(output).as_posix()
    for i,entry in enumerate(case['thrust_profile']['profiles'].values()):
        if entry.get('csv'):
            entry['csv']=archive(entry['csv'],f'thrust_{i}.csv')
    for key in ('cstar_efficiency','cf_efficiency'):
        spec=case['engine'][key]
        if spec.get('csv'):
            spec['csv']=archive(spec['csv'],key+'.csv')
    if (source_path.parent/case['files']['chamber_map']).exists():
        case['files']['chamber_map']=archive(case['files']['chamber_map'],'chamber_map.h5')
    case['torch']['standalone_purge']['output_directory']='.'
    (output/'run_config.json').write_text(json.dumps(case,indent=2)+'\n',encoding='utf-8')


def plot_results(model, rows, pressure_ladder, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import textwrap

    t=np.asarray([r['time_s'] for r in rows])
    column=lambda key:np.asarray([r[key] for r in rows])
    fig,axes=plt.subplots(4,1,figsize=(12,11),sharex=True,layout='constrained')
    for key,label in [('torch_pressure_psia','Torch chamber'),('main_chamber_backpressure_psia','Main chamber backpressure')]:
        axes[0].plot(t,column(key),label=label)
    axes[0].plot(t,pressure_ladder[:,-3],ls='--',label='FIORF outlet / fuel injector inlet')
    axes[0].set_ylabel('Pressure [psia]')
    axes[1].plot(t,1000*column('purge_mdot_kg_s'),label='Purge nitrogen into torch')
    axes[1].set_ylabel('Purge mdot [g/s]')
    axes[2].plot(t,1000*column('throat_mdot_kg_s'),label='Throat: + toward main / − backflow')
    axes[2].axhline(0,color='gray',lw=.7)
    axes[2].set_ylabel('Throat mdot [g/s]')
    for name in model.valves:
        axes[3].plot(t,column(name+'_opening'),label=name+' actual opening')
        axes[3].plot(t,column(name+'_command'),ls=':',alpha=.5,label=name+' command')
    axes[3].set_ylabel('Valve opening [0–1]')
    axes[3].set_xlabel('Simulation time [s]')
    for ax in axes:
        ax.axvline(model.release,color='gray',ls='--',lw=1)
        ax.grid(alpha=.25)
        ax.legend(fontsize=8,loc='best')
    fig.suptitle(f'Torch purge only — {model.initial_pressure/6894.757293168:g} psia imposed before '
                 f'{model.release:g} s; cold ideal gas thereafter')
    fig.savefig(output/'torch_purge_history.png',dpi=150)
    for ax in axes:
        ax.set_xlim(max(0.,model.release-.03),min(model.final_time,model.release+.2))
    fig.savefig(output/'torch_purge_transition.png',dpi=150)
    plt.close(fig)

    # Only the actual purge lifetime: omit the earlier imposed-pressure/firing phase.
    purge_valves=[p for p in model.network['purge_leg']['elements'] if p['type']=='Valve']
    starts=[when for p in purge_valves for when,opening in p['command_events_s'] if opening>0]
    purge_start=min(starts) if starts else model.release
    life=(t>=purge_start)&(t<=model.final_time)
    candidate_times={purge_start,model.final_time,(purge_start+model.final_time)/2}
    for part in purge_valves:
        for when,_ in part['command_events_s']:
            if when>=purge_start:
                stroke=max(part.get('open_time_s',0),part.get('close_time_s',0))
                candidate_times.update([when,when+stroke*.5,when+stroke,when+max(.2,stroke)])
    candidate_times.update(when for when,_ in model.main_events if when>=purge_start)
    candidate_times.update(model.settings['ladder_times_s'])
    indices=sorted({int(np.argmin(abs(t-when))) for when in candidate_times
                    if purge_start<=when<=model.final_time})
    # Keep a readable legend for long/multi-pulse cases while spanning the lifetime.
    if len(indices)>12:
        indices=[indices[i] for i in np.linspace(0,len(indices)-1,12,dtype=int)]
    fig,axes=plt.subplots(2,1,figsize=(13,10),layout='constrained')
    samples=[]
    x=np.arange(len(model.ladder_labels))
    colors=plt.get_cmap('viridis')(np.linspace(.05,.9,len(indices)))
    components=['Regulated source']+[p['name'] for p in model.parts]+[model.network['throat']['name']]
    for k,color in zip(indices,colors):
        phase='imposed Pc' if rows[k]['pressure_imposed'] else 'computed Pc'
        axes[0].plot(x,pressure_ladder[k],marker='o',color=color,label=f'{t[k]:g} s')
        if pressure_ladder[k,-2] >= float(np.min(column('main_chamber_backpressure_psia')[life]))-.05:
            axes[1].plot([0,1],pressure_ladder[k,-2:],marker='o',color=color,label=f'{t[k]:g} s')
        samples.extend({'time_s':float(t[k]),'station':label,'upstream_component':component,
                        'pressure_psia':float(p),'phase':phase}
                       for label,component,p in zip(model.ladder_labels,components,pressure_ladder[k]))
    station_titles=[model.ladder_labels[0]]+[
        (f'{part["name"]}\nOutlet' if node == part['name']+' Outlet' else f'{part["name"]} outlet\n{node}')
        for part,node in zip(model.parts,model.ladder_labels[1:-1])
    ]+[f'{model.network["throat"]["name"]} outlet\n{model.ladder_labels[-1]}']
    axes[0].set_xticks(x,['\n'.join('\n'.join(textwrap.wrap(line,22)) for line in label.splitlines())
                         for label in station_titles],fontsize=9)
    axes[0].set_title(f'Purge-only pressure ladder over its lifetime: {purge_start:g} to {model.final_time:g} s')
    axes[0].legend(fontsize=8,ncol=3,loc='best',title='Simulation time')
    # Resolve the small outward-flow pressure drop without the supply-scale axis.
    main=column('main_chamber_backpressure_psia')[life]
    positive=column('torch_over_main_pressure_psi')[life]
    padding=max(.1,float(np.max(np.maximum(positive,0)))*.25)
    axes[1].set_ylim(float(np.min(main))-padding,float(np.max(main))+max(.2,float(np.max(positive)))+padding)
    axes[1].set_xticks([0,1],[model.network['chamber']['name'],model.network['throat']['sink']])
    axes[1].set_title('Expanded torch-to-main ladder; initial low-pressure transient is outside this zoom')
    for ax in axes:
        ax.set_ylabel('Pressure [psia]')
        ax.grid(axis='y',alpha=.3)
    fig.savefig(output/'torch_pressure_ladder.png',dpi=150)
    plt.close(fig)
    write_csv(output/'torch_pressure_ladder.csv',samples)

    fig,axes=plt.subplots(3,1,figsize=(13,10),sharex=True,layout='constrained')
    for i,label in enumerate(model.ladder_labels):
        axes[0].plot(t[life],pressure_ladder[life,i],label=label)
    axes[0].set_ylabel('Station pressure [psia]')
    axes[0].legend(fontsize=8,ncol=2)
    axes[0].set_title('Purge path station pressures throughout the purge event')
    axes[1].plot(t[life],column('torch_over_main_pressure_psi')[life],color='#16776a')
    axes[1].axhline(0,color='gray',lw=.8)
    axes[1].set_ylim(-padding,max(.2,float(np.max(positive)))+padding)
    axes[1].set_ylabel('Torch above main [psi]')
    axes[1].set_title('Expanded pressure margin; negative startup transient is clipped here')
    axes[2].plot(t[life],1000*column('purge_mdot_kg_s')[life],color='#237b9b')
    axes[2].set_ylabel('Purge mdot [g/s]')
    axes[2].set_xlabel('Simulation time [s]')
    for ax in axes:
        ax.set_xlim(purge_start,model.final_time)
        ax.grid(alpha=.25)
    fig.savefig(output/'torch_purge_station_history.png',dpi=150)
    plt.close(fig)

    fig,axes=plt.subplots(2,1,figsize=(12,7),sharex=True,layout='constrained')
    margin=np.where(column('pressure_imposed')==0,column('torch_over_main_pressure_psi'),np.nan)
    for ax in axes:
        ax.plot(t,margin,color='#16776a',label='Torch pressure minus main backpressure')
        ax.axhline(0,color='gray',lw=.8)
        ax.axvline(model.release,color='gray',lw=.8,ls='--')
        ax.set_ylabel('Torch over main [psi]')
        ax.grid(alpha=.25)
    axes[0].set_title('Pressure margin after release: negative means main-to-torch backflow')
    scale=max((p['torch_over_main_pressure_psi'] for p in
               [model.fully_open_equilibrium(model.main_pressure(t)) for t in (model.release,model.final_time)]),default=1.)
    axes[1].set_ylim(-max(.02,scale*.1),max(.1,scale*1.2))
    axes[1].set_title('Same trace, expanded near zero to reveal the forward-purge pressure margin')
    axes[1].set_xlabel('Simulation time [s]')
    fig.savefig(output/'torch_pressure_margin.png',dpi=150)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description='Independent torch cold-nitrogen purge transient; no full-system run.')
    parser.add_argument('--config',help='JSON config (or HOPPER_CONFIG/default hopper_config.json).')
    parser.add_argument('--out',help='Output root; a unique UTC run directory is created inside it.')
    parser.add_argument('--no-plots',action='store_true',help='Export numeric results without rendering PNGs.')
    parser.add_argument('--check-inputs',action='store_true',help='Validate/resolve inputs without solving or writing files.')
    parser.add_argument('--diagram-only',action='store_true',help='Archive inputs and draw the full torch network without solving.')
    args=parser.parse_args()
    config,values,source_path=load_config(args.config)
    model=TorchPurgeModel(config,values)
    print(f'Purge source: {model.network["purge_leg"]["source"]}, {model.supply_pressure/6894.757293168:g} psia, {model.supply_temperature:g} K')
    print(f'Torch pressure imposed at {model.initial_pressure/6894.757293168:g} psia until {model.release:g} s.')
    print('Then solve chamber mass/energy with signed throat flow. This is cold gas, not torch combustion.')
    print('Sequence: torch.standalone_purge.sequence (hardware travel times remain in torch.network).')
    for row in model.sequence_rows():
        commands=', '.join(f'{name}={row[name]:g}' for name in model.valves)
        print(f'  {row["time_s"]:g} s: {commands}; main {row["main_chamber_pressure_psia"]:g} psia; '
              f'{row["torch_pressure_mode"]}'+(' [END]' if row['end_of_run'] else ''))
    if args.check_inputs:
        print('Inputs valid.')
        return
    root=Path(args.out).resolve() if args.out else (source_path.parent/model.settings['output_directory']).resolve()
    output=root/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    output.mkdir(parents=True,exist_ok=False)
    archive_config(config,source_path,output)
    write_csv(output/'torch_sequence.csv',model.sequence_rows())
    started=time.perf_counter()
    if not args.no_plots or args.diagram_only:
        from Backend.Hopper_Torch_Diagram import torch_diagram
        torch_diagram(model,output)
    if args.diagram_only:
        print(f'Diagram and inputs: {output}')
        return
    rows,ladder,summary=model.solve()
    write_csv(output/'torch_purge_history.csv',rows)
    ladder_rows=[{'time_s':r['time_s'],**dict(zip(model.ladder_labels,map(float,pressures)))}
                 for r,pressures in zip(rows,ladder)]
    write_csv(output/'torch_pressure_stations.csv',ladder_rows)
    write_csv(output/'torch_fully_open_equilibria.csv',summary['fully_open_purge_equilibria'])
    if not args.no_plots:
        plot_results(model,rows,ladder,output)
    summary['elapsed_s']=time.perf_counter()-started
    (output/'run_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (output/'resolved_inputs.json').write_text(json.dumps({'config':config,'gas_constant_j_kg_k':model.gas_constant,
        'gamma':model.gamma,'purge_supply_pressure_psia':model.supply_pressure/6894.757293168,
        'resolved_sequence':model.sequence_rows(),'resolved_standalone_settings':model.settings,
        'resolved_valve_commands':{name:p['command_events_s'] for name,p in model.valves.items()},
        'throat_model':'geometric area; ideal gas; bidirectional; no extra Cd','assumptions':summary['assumptions']},indent=2)+'\n',encoding='utf-8')
    print(f'Peak purge mdot: {1000*summary["peak_purge_mdot_kg_s"]:.4g} g/s')
    print(f'Peak computed torch pressure: {summary["peak_released_torch_pressure_psia"]:.5g} psia')
    print(f'Main-to-torch backflow: {1000*summary["main_chamber_backflow_kg"]:.5g} g (positive throat flow is torch-to-main).')
    for state in summary['fully_open_purge_equilibria']:
        print(f'Fully open steady purge at main {state["main_chamber_pressure_psia"]:g} psia: '
              f'torch {state["torch_pressure_psia"]:.5f} psia, margin '
              f'{state["torch_over_main_pressure_psi"]:.5f} psi, {1000*state["purge_mdot_kg_s"]:.4g} g/s; '
              f'throat {"choked" if state["throat_choked"] else "unchoked"}.')
    print(f'Finished in {summary["elapsed_s"]:.2f} s. Results: {output}')


if __name__=='__main__':
    main()
