"""Compatibility adapter. Edit hopper_config.json, not this module.

Set HOPPER_CONFIG or run Hopper.py --config path/to/config.json for another case.
All public solver constants below are loaded/derived from that JSON in SI units.
"""
from scipy.constants import *
from Hopper_Elements import (BangBangValve, CavitatingVenturi, CheckValve, COPVSource,
    Custom, GasLine, Injector, Line, Orifice, Regen, Tank, ThrottleValve, Valve,
    tube_area, tube_id)
from Hopper_Inputs import LBF
from Hopper_Config import load_config

CONFIG, _values, CONFIG_PATH = load_config()
CONFIG_SOURCE_TEXT = CONFIG_PATH.read_text(encoding='utf-8-sig')
globals().update(_values)
