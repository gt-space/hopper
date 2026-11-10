"""Torch connection drawing generated from the case's actual network definition."""
import csv
import math
from pathlib import Path
import textwrap


def torch_connections(network):
    """Component edges, retaining the shared fuel tee as an explicit graph node."""
    rows = []
    tee = network['purge_leg']['sink']
    fuel = network['fuel_leg']['elements']
    tee_index = next(i for i,p in enumerate(fuel) if p['outlet_node'] == tee)
    for key in ('fuel_leg','purge_leg','oxidizer_leg'):
        leg = network[key]
        previous = leg['source']
        for i,part in enumerate(leg['elements']):
            solved = key == 'purge_leg' or (key == 'fuel_leg' and i > tee_index)
            rows.append({'leg':key,'from':previous,'to':part['name'],'solved_in_purge':solved})
            previous = part['name']
            if key == 'fuel_leg' and i == tee_index:
                rows.append({'leg':key,'from':previous,'to':tee,'solved_in_purge':False})
                previous = tee
        rows.append({'leg':key,'from':previous,'to':leg['sink'],
                     'solved_in_purge':key != 'oxidizer_leg'})
    for start,end in ((network['chamber']['name'],network['throat']['name']),
                      (network['throat']['name'],network['throat']['sink'])):
        rows.append({'leg':'shared_outlet','from':start,'to':end,'solved_in_purge':True})
    return rows


def torch_diagram(model, output_directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon, Circle

    output = Path(output_directory)
    output.mkdir(parents=True,exist_ok=True)
    net = model.network
    fuel, purge, ox = (net[k]['elements'] for k in ('fuel_leg','purge_leg','oxidizer_leg'))
    tee = net['purge_leg']['sink']
    split = next(i for i,p in enumerate(fuel) if p['outlet_node']==tee)+1
    left_columns = max(split,len(purge))
    columns = max(left_columns+1+len(fuel)-split,len(ox))
    pitch, source_x = 2.6, 1.25
    first_x = source_x+pitch
    last_x = first_x+pitch*(columns-1)
    tee_x = first_x+pitch*left_columns
    chamber_x, throat_x, main_x = last_x+3.1,last_x+5.7,last_x+8.3
    y_fuel,y_purge,y_ox,y_chamber = 5.8,8.9,2.7,4.25
    width,height = 2.12,1.78
    ink,muted = '#22344a','#607185'
    colors={'fuel':'#b56b32','purge':'#16776a','oxidizer':'#237b9b'}
    fig,ax=plt.subplots(figsize=(22,10),facecolor='#f7f9fc')
    ax.set(xlim=(-.15,main_x+1.5),ylim=(-1.4,11.6))
    ax.axis('off')

    def arrow(points,color,active=True,bidirectional=False):
        style='-' if active else (0,(4,3))
        if len(points)>2:
            ax.plot(*zip(*points[:-1]),color=color,lw=1.7,ls=style)
        ax.add_patch(FancyArrowPatch(points[-2],points[-1],arrowstyle='<->' if bidirectional else '-|>',
                                     mutation_scale=12,lw=1.7,color=color,linestyle=style))

    def card(x,y,name,details,kind,color,active=True,tall=False):
        h=2.8 if tall else height
        ax.add_patch(FancyBboxPatch((x-width/2,y-h/2),width,h,
                     boxstyle='round,pad=0.04,rounding_size=0.1',fc='white',ec=color,
                     lw=1.4,linestyle='-' if active else (0,(4,3))))
        sx,sy=x-width/2+.2,y+h/2-.3
        if kind in ('Valve','CheckValve'):
            ax.add_patch(Polygon([(sx-.1,sy-.09),(sx+.1,sy+.09),(sx+.1,sy-.09),(sx-.1,sy+.09)],fc=color,ec=color))
        elif kind=='HeatedLine':
            for dy in (-.045,.045):
                ax.plot([sx-.11,sx+.11],[sy+dy,sy+dy],color=color,lw=1.5)
        else:
            ax.add_patch(Circle((sx,sy),.09,fc='white',ec=color,lw=1.4))
        ax.text(x-width/2+.39,sy,'\n'.join(textwrap.wrap(name,20)),va='center',
                fontsize=9,weight='bold',color=ink)
        ax.text(x-width/2+.12,y+h/2-.65,'\n'.join(details),va='top',fontsize=8.5,
                color=muted,linespacing=1.55)

    def details(part):
        if part['type']=='HeatedLine':
            od,wall,diameter=(part.get(k) for k in ('outer_diameter_in','wall_thickness_in','inner_diameter_in'))
            if diameter is None and od is not None and wall is not None:
                diameter=od-2*wall
            fmt=lambda v:'unset' if v is None else f'{v:g}'
            return [f'L {fmt(part.get("length_in"))} in',f'OD {fmt(od)} / wall {fmt(wall)} in',
                    f'ID {fmt(diameter)} in','Input dT vs torch time']
        area=part.get('cda_mm2')
        lines=[f'CdA {"unset" if area is None else f"{area:g}"} mm²']
        if part['type']=='Valve':
            lines += [f'Open {part.get("open_time_s",0):g} s',f'Close {part.get("close_time_s",0):g} s']
        elif part['type']=='CheckValve':
            lines += ['Forward flow only','Cracking dP not modelled']
        else:
            lines += [part['type']]
        return lines

    ax.text(0,11.1,'TORCH FLUIDS  /  CONFIGURED CONNECTIONS',fontsize=19,weight='bold',color=ink)
    ax.text(0,10.62,'Solid: solved cold-gas purge path     |     Dashed: configured propellant feeds, not solved in purge-only mode',
            fontsize=10,color=muted)
    for key,label,y,color_key in (('purge_leg','REGULATED GN2 PURGE',y_purge,'purge'),
                                 ('fuel_leg','FUEL FEED',y_fuel,'fuel'),
                                 ('oxidizer_leg','OXIDIZER FEED',y_ox,'oxidizer')):
        leg,color=net[key],colors[color_key]
        active=key=='purge_leg'
        ax.text(0,y+1.16,label,fontsize=10,weight='bold',color=color)
        source_details=([f'{model.supply_pressure/6894.757293168:.5g} psia',f'{model.supply_temperature:g} K',
                         'Shared main-purge supply'] if active else ['Upstream of main','throttle valve','Boundary; feed not solved'])
        card(source_x,y,leg['source'],source_details,'Source',color,active)

    fuel_x = ([first_x+pitch*i for i in range(split)]
              +[tee_x+pitch*(i+1) for i in range(len(fuel)-split)])
    purge_x=[first_x+pitch*i for i in range(len(purge))]
    ox_x=[first_x+pitch*i for i in range(len(ox))]
    for parts,positions,y,color,key in ((fuel,fuel_x,y_fuel,colors['fuel'],'fuel'),
                                       (purge,purge_x,y_purge,colors['purge'],'purge'),
                                       (ox,ox_x,y_ox,colors['oxidizer'],'ox')):
        previous=source_x
        for i,(part,x) in enumerate(zip(parts,positions)):
            active=key=='purge' or (key=='fuel' and i>=split)
            line_color=colors['purge'] if active else color
            card(x,y,part['name'],details(part),part['type'],line_color,active)
            if key=='fuel' and i==split:
                arrow([(previous+width/2,y),(tee_x,y)],color,False)
                arrow([(tee_x,y),(x-width/2,y)],colors['purge'])
            else:
                arrow([(previous+width/2,y),(x-width/2,y)],line_color,active)
            previous=x
    arrow([(purge_x[-1]+width/2,y_purge),(tee_x,y_purge),(tee_x,y_fuel)],colors['purge'])
    ax.add_patch(Circle((tee_x,y_fuel),.075,fc=colors['purge'],ec='white',zorder=6))
    ax.text(tee_x,y_fuel-.28,f'{tee}\nShared purge tee',ha='center',va='top',fontsize=9,color=ink)

    for x,y,port_y,color,active in ((fuel_x[-1],y_fuel,y_chamber+.55,colors['purge'],True),
                                   (ox_x[-1],y_ox,y_chamber-.55,colors['oxidizer'],False)):
        rail=chamber_x-width/2-.5
        arrow([(x+width/2,y),(rail,y),(rail,port_y),(chamber_x-width/2,port_y)],color,active)
    chamber=net['chamber']
    card(chamber_x,y_chamber,chamber['name'],[f'Volume {chamber["volume_l"]:g} L',
         f'Initial {chamber["initial_pressure_psia"]:g} psia',f'Initial {chamber["initial_temperature_k"]:g} K',
         f'Imposed until {model.release:g} s','Then mass + energy','No combustion'], 'Chamber',ink,tall=True)
    throat=net['throat']
    area=throat['area_in2']
    card(throat_x,y_chamber,throat['name'],[f'Area {area:g} in²',f'Diameter {math.sqrt(4*area/math.pi):.5g} in',
         'Geometric area; no Cd','Bidirectional gas flow'],'Throat',ink)
    pressures=' → '.join(f'{p:g}' for _,p in model.main_events)
    card(main_x,y_chamber,throat['sink'],['Backpressure boundary',pressures+' psia',
         f'{model.settings["main_chamber_gas_temperature_k"]:g} K cold-gas assumption'],'Source',ink)
    arrow([(chamber_x+width/2,y_chamber),(throat_x-width/2,y_chamber)],ink,bidirectional=True)
    arrow([(throat_x+width/2,y_chamber),(main_x-width/2,y_chamber)],ink,bidirectional=True)
    ax.text((throat_x+main_x)/2,y_chamber-1.35,'Reverse flow allowed\nif main pressure > torch',ha='center',fontsize=9,color=muted)
    ax.text(0,.65,'VERIFY THE TEE:  purge joins after FIORF and before the SAME fuel injector used by the fuel feed.',
            fontsize=11,weight='bold',color=ink)
    ax.text(0,.13,'Throat connects torch to main chamber. Main backpressure is imposed; torch pressure is solved after release.',
            fontsize=10,color=muted)
    ax.text(0,-.42,'Timing: torch.standalone_purge.sequence     |     Hardware: torch.network     |     No liquid displacement or hot-gas chemistry.',
            fontsize=10,color=muted)
    fig.tight_layout(pad=.7)
    for extension in ('svg','png'):
        fig.savefig(output/f'torch_component_diagram.{extension}',dpi=160,facecolor=fig.get_facecolor())
    plt.close(fig)
    with (output/'torch_connections.csv').open('w',newline='',encoding='utf-8') as file:
        rows=torch_connections(net)
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output/'torch_component_diagram.svg'
