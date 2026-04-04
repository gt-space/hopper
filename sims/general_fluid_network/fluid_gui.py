import sys
import inspect
import json  # <--- NEW
from PyQt5.QtWidgets import (QApplication, QMainWindow, QGraphicsScene, QGraphicsView, 
                             QDockWidget, QListWidget, QVBoxLayout, QWidget, QLabel, 
                             QFormLayout, QLineEdit, QPushButton, QGraphicsItem, 
                             QGraphicsLineItem, QInputDialog, QMessageBox, QGraphicsPathItem, QFileDialog, QTableWidgetItem, QComboBox,
                             QTableWidget, QHeaderView, QGraphicsDropShadowEffect, QDialog)
from PyQt5.QtCore import Qt, QMimeData, QPointF, QRectF, QLineF
from PyQt5.QtGui import QDrag, QPainter, QPen, QBrush, QColor, QPainterPath, QImage, QPainterPathStroker

# --- IMPORT YOUR SIMULATION LIBRARY ---
import general_fluid_network as gfn
import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
#  CONFIGURATION & EXTENSIBILITY
# =============================================================================
# Define available components here.
# Keys: Class names from your library.
# Values: Dictionary of default __init__ parameters.
# TO ADD A NEW OBJECT: Just add it to this dictionary!
COMPONENT_DEFS = {
    "Nodes": {
        "Node": {
            "fluid": "Water", "m": 1.0, "V": 1.0, "T": 293.15, 
            "name": "node", "type": "m"
        },
        "Tank": {
            "V_total_L": 100.0, "fluid_liq": "Water", "m_liq": 10.0, "T_liq": 293.15,
            "fluid_ullage": "Nitrogen", "P_ullage": 101325.0, "T_ullage": 293.15,
            "radius": 0.5, "htc": 50.0, "name": "tank"
        },
        "Ambient": {
            "fluid": "Air", "P": 101325.0, "T": 293.15, "name": "ambient"
        },
        "Manifold": {"fluid": "Nitrogen", "m": 0.1, "V": 1.0, "T": 293.15, "name": "manifold"},
        "Engine": {
            "fuel": "Kerosene",
            "oxidizer": "LOX",
            "mdot_ox": 2.6,
            "mdot_fuel": 1.0,
            "Pc": 3e6,
            "eta_cstar": 0.95,
            "At": 0.002026,
            "Ae": 0.01013,
            "Pa": 101325.0,
            "name": "engine"
        },
    },
    "Connections": {
        "Connection": {
            "CdA": 0.0001, "qdot": 0.0, "location": 0.0, 
            "normal_state": 1, "checking": 1, "name": "pipe"
        },
        "Valve": {
            "CdA": 0.0001, "qdot": 0.0, "location": 0.0, 
            "normal_state": 0, "checking": 1, "name": "valve"
        },
        "CheckValve": {
            "CdA": 0.0001, "qdot": 0.0, "location": 0.0, 
            "normal_state": 1, "checking": 1, "name": "check_valve"
        },
        "Regulator": {
            "CdA": 0.0001, "set_pressure": 200000.0, "qdot": 0.0, 
            "location": 0.0, "normal_state": 1, "name": "regulator"
        },
        "ThrottleValve": {
            "CdA_max": 0.0001, "qdot": 0.0, 
            "location": 0.0, "normal_state": 0, "checking": 1, 
            "name": "throttle"
        },
    
        "BangBang": {
            "CdA": 0.0001, "target_pressure": 150000.0, "hysteresis": 5000.0,
            "qdot": 0.0, "location": 0.0, "normal_state": 1, 
            "checking": 1, "name": "bang_bang"
        },
        "Line": {
            "ID": 0.001, "length": 0.1, "roughness": 0.0002, "qdot": 0.0,
            "location": 0.0, "normal_state": 1, "checking": 1, "name": "line"
        },
        "Series": {
            "connections": [], "name": "series"
        },
        "Engine": {
            "fuel": "n-Dodecane",
            "oxidizer": "Oxygen",
            "eta_cstar": 0.92,
            "At": 0.002,
            "Ae": 0.01,
            "Pa": 101325.0,
            "name": "engine"
        }
    }
}

# =============================================================================
#  FLUID PROPERTY LOOKUP (REFPROP -> CoolProp)
# =============================================================================

def _normalize_prop_key(key):
    if not isinstance(key, str):
        return ""
    return key.strip().upper()

def _try_props_refprop_then_coolprop(fluid, out_key, in1_key, in1_val, in2_key, in2_val):
    """
    Try REFPROP backend first (via CoolProp), then fallback to CoolProp default backend.
    Units: SI (T K, P Pa, D kg/m3, H J/kg, S J/kg/K, Q 0-1).
    """
    try:
        from CoolProp.CoolProp import PropsSI
    except Exception as e:
        raise RuntimeError(f"CoolProp not available: {e}")

    out_key = _normalize_prop_key(out_key)
    in1_key = _normalize_prop_key(in1_key)
    in2_key = _normalize_prop_key(in2_key)

    if not fluid or not out_key or not in1_key or not in2_key:
        raise ValueError("Fluid and property keys are required.")

    # Try REFPROP backend first through CoolProp
    refprop_fluid = f"REFPROP::{fluid}"
    try:
        return PropsSI(out_key, in1_key, in1_val, in2_key, in2_val, refprop_fluid)
    except Exception:
        # Fallback to CoolProp default backend
        return PropsSI(out_key, in1_key, in1_val, in2_key, in2_val, fluid)

# =============================================================================
#  UNITS MAP (UI LABELS)
# =============================================================================

UNITS_MAP = {
    # Common state variables
    "T": "K",
    "P": "Pa",
    "D": "kg/m3",
    "H": "J/kg",
    "S": "J/kg/K",
    "Q": "-",
    "m": "kg",
    "V": "m3",
    # Component params
    "V_total_L": "L",
    "m_liq": "kg",
    "T_liq": "K",
    "P_ullage": "Pa",
    "T_ullage": "K",
    "radius": "m",
    "htc": "W/m2/K",
    "CdA": "m2",
    "CdA_max": "m2",
    "qdot": "W",
    "location": "m",
    "set_pressure": "Pa",
    "hysteresis": "Pa",
    #Engine units
    "Pc": "Pa",
    "Pa": "Pa",
    "mdot_ox": "kg/s",
    "mdot_fuel": "kg/s",
    "At": "m²",
    "Ae": "m²",
}

def _label_with_units(key):
    unit = UNITS_MAP.get(key)
    return f"{key} [{unit}]" if unit else key

# =============================================================================
#  UI CONSTANTS
# =============================================================================

GLOW_COLOR = QColor(0, 200, 255, 180)
COMMON_FLUIDS = ["Water", "Air", "Nitrogen", "Oxygen", "Methane", "Propane", "CO2"]

# =============================================================================
#  GRAPHICS ITEMS (The Visuals)
# =============================================================================

class NodeItem(QGraphicsItem):
    def __init__(self, node_type, params, pos):
        super().__init__()
        self.node_type = node_type
        self.params = params.copy()
        self.setPos(pos)
        self.setFlags(QGraphicsItem.ItemIsMovable | QGraphicsItem.ItemIsSelectable | QGraphicsItem.ItemSendsGeometryChanges)
        self.sim_instance = None 
        
        self.width = 80
        self.height = 60
        self.color = QColor("#3498db") if "Tank" not in node_type else QColor("#e67e22")
        self._glow = QGraphicsDropShadowEffect()
        self._glow.setBlurRadius(25)
        self._glow.setColor(GLOW_COLOR)
        self._glow.setOffset(0, 0)
        self._glow.setEnabled(False)
        self.setGraphicsEffect(self._glow)

    def boundingRect(self):
        return QRectF(-self.width/2, -self.height/2, self.width, self.height)

    def paint(self, painter, option, widget):
        painter.setPen(QPen(Qt.black, 2 if self.isSelected() else 1))
        painter.setBrush(QBrush(self.color))
        painter.drawRoundedRect(self.boundingRect(), 5, 5)
        
        painter.setPen(Qt.white)
        name = self.params.get("name", self.node_type)
        painter.drawText(self.boundingRect(), Qt.AlignCenter, name)

    def itemChange(self, change, value):
        if change == QGraphicsItem.ItemSelectedHasChanged:
            if self._glow:
                self._glow.setEnabled(bool(self.isSelected()))
        if change == QGraphicsItem.ItemPositionChange:
            scene = self.scene()
            if scene:
                scene.update_connections(self)
        return super().itemChange(change, value)

    # --- THE FIX IS HERE ---
    def mousePressEvent(self, event):
        # If CTRL is pressed, IGNORE this event.
        # This lets the event "bubble up" to the View, which handles the connection drawing.
        if event.modifiers() == Qt.ControlModifier:
            event.ignore()
        else:
            # Otherwise, let the Node handle it (Selection / Moving)
            super().mousePressEvent(event)

class ConnectionLine(QGraphicsPathItem):
    def __init__(self, start_item, end_item, conn_type, params):
        super().__init__()
        self.start_item = start_item
        self.end_item = end_item
        self.conn_type = conn_type
        self.params = params.copy()
        self.sim_instance = None
        
        self.setFlags(QGraphicsItem.ItemIsSelectable)
        self.setZValue(-1) 
        self.update_position()

    def update_position(self):
        line = QLineF(self.start_item.pos(), self.end_item.pos())
        path = QPainterPath()
        path.moveTo(line.p1())
        path.lineTo(line.p2())
        self.setPath(path)
        
        # Standard Visual Pen
        pen = QPen(Qt.darkGray, 3)
        if self.isSelected():
            pen.setColor(Qt.red)
            pen.setStyle(Qt.DashLine)
        self.setPen(pen)

    def shape(self):
        # --- THE FIX FOR CLICKING ---
        # Create a "Stroker" that creates a shape from the path
        # We make this shape much thicker (20px) than the visible line
        path = self.path()
        path_stroker = QPainterPathStroker()
        path_stroker.setWidth(20) # 20px hit area
        return path_stroker.createStroke(path)

    def paint(self, painter, option, widget):
        if self.isSelected():
            glow_pen = QPen(GLOW_COLOR, 10, Qt.SolidLine, Qt.RoundCap)
            painter.setPen(glow_pen)
            painter.setBrush(Qt.NoBrush)
            painter.drawPath(self.path())
        super().paint(painter, option, widget)
        
        # Draw label box
        mid_point = self.path().pointAtPercent(0.5)
        painter.setBrush(Qt.white)
        painter.setPen(Qt.black)
        # Center the rect
        rect = QRectF(mid_point.x()-15, mid_point.y()-10, 30, 20)
        painter.drawRect(rect)
        
        painter.setFont(painter.font())
        label = self.conn_type[:3].upper()
        painter.drawText(rect, Qt.AlignCenter, label)
# =============================================================================
#  THE EDITOR CANVAS
# =============================================================================

class FluidScene(QGraphicsScene):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.line_preview = None
        self.connect_start_item = None

    def update_connections(self, node_item):
        for item in self.items():
            if isinstance(item, ConnectionLine):
                if item.start_item == node_item or item.end_item == node_item:
                    item.update_position()

class FluidEditor(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("VibeFlow - Fluid Network Designer")
        self.resize(1200, 800)

        # 1. Main Canvas
        self.scene = FluidScene()
        self.view = QGraphicsView(self.scene)
        self.view.setRenderHint(QPainter.Antialiasing)
        self.setCentralWidget(self.view)

        # Property calc state
        self._prop_calc_widgets = {}
        self._prop_calc_last = {
            "p1": "T",
            "p2": "P",
            "out": "D",
        }
        self._prop_calc_item = None

        # 2. Component Library (Sidebar)
        self.create_library_dock()
        
        # 3. Property Editor (Sidebar)
        self.create_property_dock()

        # 4. Fluid Property Calc (Bottom Left)
        self.create_prop_calc_dock()

        # 5. Simulation Controls (Top)
        self.create_sim_dock()

        # Event wiring
        self.scene.selectionChanged.connect(self.on_selection_changed)
        self.view.setAcceptDrops(True)

    def create_library_dock(self):
        dock = QDockWidget("Library", self)
        dock.setObjectName("LibraryDock")
        list_widget = QListWidget()
        
        # Populate from Config
        for category, items in COMPONENT_DEFS.items():
            # Add header? maybe later. Just flat list for now or colored
            for name in items.keys():
                list_widget.addItem(name)
        
        list_widget.setDragEnabled(True)
        list_widget.startDrag = self.start_drag_library # Custom drag handler
        dock.setWidget(list_widget)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)
        self.library_list = list_widget
        self.library_dock = dock

    def start_drag_library(self, supportedActions):
        item = self.library_list.currentItem()
        drag = QDrag(self)
        mime = QMimeData()
        mime.setText(item.text())
        drag.setMimeData(mime)
        drag.exec_(Qt.CopyAction)

    def create_property_dock(self):
        dock = QDockWidget("Properties", self)
        self.prop_widget = QWidget()
        self.prop_layout = QFormLayout()
        self.prop_widget.setLayout(self.prop_layout)
        dock.setWidget(self.prop_widget)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)
        self.current_editing_item = None

    def create_prop_calc_dock(self):
        dock = QDockWidget("Fluid Property Calc", self)
        widget = QWidget()
        layout = QFormLayout()

        prop_keys = ["T", "P", "D", "H", "S", "Q"]

        self.calc_cmb_fluid = QComboBox()
        self.calc_cmb_fluid.setEditable(True)
        self.calc_cmb_fluid.addItems(COMMON_FLUIDS)
        self.calc_cmb_fluid.setCurrentText("Water")
        layout.addRow("fluid", self.calc_cmb_fluid)

        self.calc_cmb_p1 = QComboBox()
        self.calc_cmb_p1.addItems(prop_keys)
        self.calc_cmb_p1.setCurrentText(self._prop_calc_last["p1"])
        self.calc_inp_p1 = QLineEdit("")
        self.calc_lbl_p1v = QLabel(_label_with_units(self.calc_cmb_p1.currentText()))
        layout.addRow("prop1", self.calc_cmb_p1)
        layout.addRow(self.calc_lbl_p1v, self.calc_inp_p1)

        self.calc_cmb_p2 = QComboBox()
        self.calc_cmb_p2.addItems(prop_keys)
        self.calc_cmb_p2.setCurrentText(self._prop_calc_last["p2"])
        self.calc_inp_p2 = QLineEdit("")
        self.calc_lbl_p2v = QLabel(_label_with_units(self.calc_cmb_p2.currentText()))
        layout.addRow("prop2", self.calc_cmb_p2)
        layout.addRow(self.calc_lbl_p2v, self.calc_inp_p2)

        self.calc_cmb_out = QComboBox()
        self.calc_cmb_out.addItems(prop_keys)
        self.calc_cmb_out.setCurrentText(self._prop_calc_last["out"])
        self.calc_out_val = QLineEdit("")
        self.calc_out_val.setReadOnly(True)
        self.calc_lbl_outv = QLabel(_label_with_units(self.calc_cmb_out.currentText()))
        layout.addRow("target", self.calc_cmb_out)
        layout.addRow(self.calc_lbl_outv, self.calc_out_val)

        btn_calc = QPushButton("Compute")
        btn_calc.clicked.connect(self.compute_fluid_property)
        layout.addRow("", btn_calc)

        widget.setLayout(layout)
        dock.setWidget(widget)
        if hasattr(self, "library_dock") and self.library_dock:
            self.addDockWidget(Qt.LeftDockWidgetArea, dock)
            self.splitDockWidget(self.library_dock, dock, Qt.Vertical)
        else:
            self.addDockWidget(Qt.LeftDockWidgetArea, dock)

        self._prop_calc_widgets = {
            "fluid": self.calc_cmb_fluid,
            "p1": self.calc_cmb_p1,
            "p1_val": self.calc_inp_p1,
            "p2": self.calc_cmb_p2,
            "p2_val": self.calc_inp_p2,
            "out": self.calc_cmb_out,
            "out_val": self.calc_out_val,
            "lbl_p1_val": self.calc_lbl_p1v,
            "lbl_p2_val": self.calc_lbl_p2v,
            "lbl_out_val": self.calc_lbl_outv,
        }

        self.calc_cmb_p1.currentTextChanged.connect(self.update_prop_calc_labels)
        self.calc_cmb_p2.currentTextChanged.connect(self.update_prop_calc_labels)
        self.calc_cmb_out.currentTextChanged.connect(self.update_prop_calc_labels)

    def create_sim_dock(self):
        dock = QDockWidget("Simulation", self)
        widget = QWidget()
        layout = QVBoxLayout()
        
        self.lbl_time = QLabel("Duration (s):")
        self.inp_time = QLineEdit("5.0")
        self.lbl_dt = QLabel("Timestep (s):")
        self.inp_dt = QLineEdit("0.01")
        
        btn_run = QPushButton("RUN SIMULATION")
        btn_run.setStyleSheet("background-color: #2ecc71; color: white; font-weight: bold; padding: 10px;")
        btn_run.clicked.connect(self.run_simulation)
        self.btn_run = btn_run
        
        layout.addWidget(self.lbl_time)
        layout.addWidget(self.inp_time)
        layout.addWidget(self.lbl_dt)
        layout.addWidget(self.inp_dt)
        layout.addWidget(btn_run)
        layout.addStretch()
        
        widget.setLayout(layout)
        dock.setWidget(widget)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)

    def _show_sim_running_dialog(self):
        dlg = QDialog(self)
        dlg.setWindowTitle("Simulation Running")
        dlg.setModal(False)
        dlg.setWindowFlag(Qt.WindowContextHelpButtonHint, False)
        layout = QVBoxLayout()
        layout.addWidget(QLabel("Simulation running..."))
        dlg.setLayout(layout)
        dlg.resize(240, 80)
        dlg.show()
        QApplication.processEvents()
        return dlg

    # --- DRAG AND DROP ONTO CANVAS ---
    def dragEnterEvent(self, event): # On QMainWindow
        pass # Handled by view usually, but we need to override view events

    # We need to monkey-patch or subclass View for proper drop handling
    # For brevity in a single file, we will wire the view events here:
    def setup_view_events(self):
        # This is a bit hacky for a single script, standard way is subclassing QGraphicsView
        pass

    # --- CONNECTION LOGIC (CTRL + DRAG) ---
    # We'll use mouse overrides on the VIEW
    # But since we didn't subclass View, let's just do it in the Scene or check modifiers
    pass 
    
# To keep this clean, let's subclass View properly
class DesignView(QGraphicsView):
    def __init__(self, scene, parent_window):
        super().__init__(scene)
        self.parent_window = parent_window
        self.setAcceptDrops(True)
        self.temp_line = None
        self.start_node = None

    # --- DRAG & DROP EVENTS ---
    def dragEnterEvent(self, event):
        if event.mimeData().hasText(): event.accept()
        else: event.ignore()

    def dragMoveEvent(self, event):
        event.accept()

    def dropEvent(self, event):
        name = event.mimeData().text()
        pos = self.mapToScene(event.pos())
        if name in COMPONENT_DEFS["Nodes"]:
            defaults = COMPONENT_DEFS["Nodes"][name]
            item = NodeItem(name, defaults, pos)
            self.scene().addItem(item)
        else:
            QMessageBox.information(self, "Info", "Connections must be drawn between two Nodes.\nHold CTRL and drag from one Node to another.")
        event.accept()

    # --- MOUSE EVENTS ---
    def mousePressEvent(self, event):
        if event.modifiers() == Qt.ControlModifier:
            pos = self.mapToScene(event.pos())
            items = self.scene().items(pos)
            for item in items:
                if isinstance(item, NodeItem):
                    self.start_node = item
                    self.temp_line = QGraphicsLineItem(QLineF(item.pos(), item.pos()))
                    self.temp_line.setPen(QPen(Qt.black, 2, Qt.DashLine))
                    self.scene().addItem(self.temp_line)
                    break
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self.temp_line and self.start_node:
            new_pos = self.mapToScene(event.pos())
            line = QLineF(self.start_node.pos(), new_pos)
            self.temp_line.setLine(line)
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if self.temp_line:
            self.scene().removeItem(self.temp_line)
            self.temp_line = None
            pos = self.mapToScene(event.pos())
            items = self.scene().items(pos)
            end_node = None
            for item in items:
                if isinstance(item, NodeItem) and item != self.start_node:
                    end_node = item
                    break
            if end_node:
                self.parent_window.create_connection(self.start_node, end_node)
            self.start_node = None
        super().mouseReleaseEvent(event)
    
    # --- PLOTTING LOGIC (Updated to Match GFN Exactly) ---
    #chat did this so lowkey might've messed other things up
    def mouseDoubleClickEvent(self, event):
        item = self.scene().itemAt(self.mapToScene(event.pos()), self.transform())
        
        if item and item.parentItem():
            item = item.parentItem()

        if isinstance(item, (NodeItem, ConnectionLine)) and item.sim_instance:

            instance = item.sim_instance

            # --- Engine special case ---
            if instance.__class__.__name__ == "Engine":
                msg = QMessageBox()
                msg.setWindowTitle("Engine Performance")

                msg.setText(instance.get_summary_text())

                msg.exec_()
                return

            # --- All other components behave normally ---
            self.plot_detailed_results(item)

        super().mouseDoubleClickEvent(event)

    def plot_detailed_results(self, item):
        """
        Plots 6-panel graphs matching general_fluid_network.py behavior.
        Units="E" (English) enforced to match user request.
        """
        obj = item.sim_instance
        hist = obj.history
        t = hist['time']
        
        if len(t) == 0:
            QMessageBox.warning(self, "No Data", "Run simulation first!")
            return

        # Setup 2x3 Grid
        fig, axs = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
        fig.suptitle(f"Results: {item.params.get('name', 'Object')}", fontsize=14)
        axs = axs.flatten()

        # --- NODES (Matches plot_nodes_overlay) ---
        if isinstance(item, NodeItem):
            # 1. Pressure (psi)
            axs[0].plot(t, np.array(hist['P']) / 6894.75729, color='blue')
            axs[0].set_ylabel("Pressure [psi]")
            
            # 2. Temperature (F)
            axs[1].plot(t, (np.array(hist['T']) - 273.15) * 1.8 + 32, color='red')
            axs[1].set_ylabel("Temperature [F]")
            
            # 3. Mass (kg) - GFN keeps this in kg for "E" units
            axs[2].plot(t, hist['m'], color='green')
            axs[2].set_ylabel("Mass [kg]")
            
            # 4. Density (kg/m3)
            axs[3].plot(t, hist['d'], color='purple')
            axs[3].set_ylabel("Density [kg/m³]")
            
            # 5. Quality (0-1)
            axs[4].plot(t, hist['Q'], color='orange')
            axs[4].set_ylabel("Quality [-]")
            axs[4].set_ylim(0, 1.1)
            
            # 6. Fill Level (0-1)
            axs[5].plot(t, hist['fill_level'], color='teal')
            axs[5].set_ylabel("Fill Level [-]")
            axs[5].set_xlabel("Time [s]")

        # --- CONNECTIONS (Matches plot_connections_overlay) ---
        elif isinstance(item, ConnectionLine):
            # 1. Mass Flow (kg/s) - GFN keeps this in kg/s for "E" units
            axs[0].plot(t, hist['mdot'], color='black')
            axs[0].set_ylabel("mdot [kg/s]")
            
            # 2. Delta Pressure (psi)
            axs[1].plot(t, np.array(hist['dP']) / 6894.75729, color='orange')
            axs[1].set_ylabel("dP [psi]")
            
            # 3. CdA (mm^2) - GFN scales by 1e6
            axs[2].plot(t, np.array(hist['CdA']) * 1000000, color='blue')
            axs[2].set_ylabel("CdA [mm²]")
            
            # 4. Enthalpy Flow (J/s)
            axs[3].plot(t, hist['Hdot'], color='red')
            axs[3].set_ylabel("Hdot [J/s]")
            
            # 5. Quality
            # Filter None values for plotting
            Q_clean = [q if q is not None else 0 for q in hist['Q']]
            axs[4].plot(t, Q_clean, color='green')
            axs[4].set_ylabel("Quality [-]")
            axs[4].set_ylim(0, 1.1)
            
            # 6. State (0-1)
            axs[5].plot(t, hist['state'], color='gray', linestyle='--')
            axs[5].set_ylabel("State [-]")
            axs[5].set_xlabel("Time [s]")

        for ax in axs:
            ax.grid(True, alpha=0.3)
            ax.legend([item.params.get("name", "Data")])

        plt.tight_layout(rect=[0, 0, 1, 0.95]) # Adjust for suptitle
        plt.show()

class MainApp(FluidEditor):
    def __init__(self):
        super().__init__()
        self.view = DesignView(self.scene, self)
        self.view.setRenderHint(QPainter.Antialiasing)
        self.setCentralWidget(self.view)
        
        # --- 1. MENU BAR (Save/Load) ---
        menubar = self.menuBar()
        file_menu = menubar.addMenu('&File')
        file_menu.addAction('&Save Network', self.save_network, 'Ctrl+S')
        file_menu.addAction('&Load Network', self.load_network, 'Ctrl+O')

        # --- 2. INITIALIZE DOCKS ---
        self.create_actions_dock()

    def create_sim_dock(self):
        # Overriding parent to include Save/Load buttons
        dock = QDockWidget("Simulation & Files", self)
        widget = QWidget()
        layout = QVBoxLayout()
        
        self.lbl_time = QLabel("Duration (s):")
        self.inp_time = QLineEdit("5.0")
        self.lbl_dt = QLabel("Timestep (s):")
        self.inp_dt = QLineEdit("0.01")
        
        btn_run = QPushButton("RUN SIMULATION")
        btn_run.setStyleSheet("background-color: #2ecc71; color: white; font-weight: bold; padding: 10px;")
        btn_run.clicked.connect(self.run_simulation)
        self.btn_run = btn_run
        
        # File Operations
        btn_save = QPushButton("Save Network")
        btn_save.clicked.connect(self.save_network)
        btn_load = QPushButton("Load Network")
        btn_load.clicked.connect(self.load_network)
        
        layout.addWidget(self.lbl_time)
        layout.addWidget(self.inp_time)
        layout.addWidget(self.lbl_dt)
        layout.addWidget(self.inp_dt)
        layout.addWidget(btn_run)
        
        layout.addSpacing(20)
        layout.addWidget(QLabel("<b>File Operations:</b>"))
        layout.addWidget(btn_save)
        layout.addWidget(btn_load)
        layout.addStretch()
        
        widget.setLayout(layout)
        dock.setWidget(widget)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)

    def create_actions_dock(self):
        dock = QDockWidget("Script Actions", self)
        widget = QWidget()
        layout = QVBoxLayout()
        
        # Table Setup
        self.action_table = QTableWidget(0, 3)
        self.action_table.setHorizontalHeaderLabels(["Time (s)", "Component", "New State"])
        header = self.action_table.horizontalHeader()
        header.setSectionResizeMode(1, QHeaderView.Stretch)
        
        btn_add = QPushButton("Add Action")
        btn_add.clicked.connect(self.add_action_row)
        btn_remove = QPushButton("Remove Selected")
        btn_remove.clicked.connect(self.remove_action_row)
        
        layout.addWidget(self.action_table)
        layout.addWidget(btn_add)
        layout.addWidget(btn_remove)
        
        widget.setLayout(layout)
        dock.setWidget(widget)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)

    def add_action_row(self):
        row = self.action_table.rowCount()
        self.action_table.insertRow(row)
        
        # Time
        self.action_table.setItem(row, 0, QTableWidgetItem("1.0"))
        
        # Component Dropdown
        combo = QComboBox()
        conn_items = [i for i in self.scene.items() if isinstance(i, ConnectionLine)]
        for item in conn_items:
            name = item.params.get("name", "Unnamed")
            combo.addItem(name, userData=item) 
            
            # NEW: Add sub-connections to the dropdown if it's a Series!
            if item.conn_type == "Series":
                for idx, sub_conn in enumerate(item.params.get("connections", [])):
                    sub_name = sub_conn["params"].get("name", f"Sub {idx+1}")
                    # userData stores a tuple: (Parent Series Item, Index of Sub-Component)
                    combo.addItem(f"  ↳ {sub_name}", userData=(item, idx))
                    
        self.action_table.setCellWidget(row, 1, combo)
        
        # State
        self.action_table.setItem(row, 2, QTableWidgetItem("0"))

    def remove_action_row(self):
        curr_row = self.action_table.currentRow()
        if curr_row >= 0:
            self.action_table.removeRow(curr_row)

    # --- SAVE / LOAD LOGIC ---
    def save_network(self):
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Network", "", "JSON Files (*.json)")
        if not file_path: return

        data = {
            "version": 1.1,
            "settings": {
                "duration": self.inp_time.text(),
                "dt": self.inp_dt.text()
            },
            "nodes": [],
            "connections": [],
            "actions": [] 
        }

        # 1. Save Nodes
        node_items = [i for i in self.scene.items() if isinstance(i, NodeItem)]
        item_to_id = {}
        for idx, item in enumerate(node_items):
            item_to_id[item] = idx
            data["nodes"].append({
                "id": idx,
                "type": item.node_type,
                "x": item.pos().x(),
                "y": item.pos().y(),
                "params": item.params
            })

        # 2. Save Connections
        conn_items = [i for i in self.scene.items() if isinstance(i, ConnectionLine)]
        for item in conn_items:
            if item.start_item in item_to_id and item.end_item in item_to_id:
                data["connections"].append({
                    "type": item.conn_type,
                    "start_id": item_to_id[item.start_item],
                    "end_id": item_to_id[item.end_item],
                    "params": item.params
                })

        # 3. Save Actions Table
        for row in range(self.action_table.rowCount()):
            t_item = self.action_table.item(row, 0)
            state_item = self.action_table.item(row, 2)
            combo = self.action_table.cellWidget(row, 1)
            
            if t_item and state_item and combo:
                comp_name = combo.currentText()
                data["actions"].append({
                    "time": t_item.text(),
                    "component": comp_name,
                    "state": state_item.text()
                })

        try:
            with open(file_path, 'w') as f:
                json.dump(data, f, indent=4)
            QMessageBox.information(self, "Success", "Network and Scripts saved successfully!")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save: {str(e)}")

    def load_network(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Load Network", "", "JSON Files (*.json)")
        if not file_path: return

        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            self.scene.clear()
            self.action_table.setRowCount(0) 
            
            # 1. Settings
            if "settings" in data:
                self.inp_time.setText(str(data["settings"].get("duration", 5.0)))
                self.inp_dt.setText(str(data["settings"].get("dt", 0.01)))

            # 2. Nodes
            id_to_item = {}
            for n_data in data["nodes"]:
                pos = QPointF(n_data["x"], n_data["y"])
                item = NodeItem(n_data["type"], n_data["params"], pos)
                self.scene.addItem(item)
                id_to_item[n_data["id"]] = item

            # 3. Connections
            loaded_connections = [] 
            for c_data in data["connections"]:
                start = id_to_item.get(c_data["start_id"])
                end = id_to_item.get(c_data["end_id"])
                if start and end:
                    item = ConnectionLine(start, end, c_data["type"], c_data["params"])
                    self.scene.addItem(item)
                    loaded_connections.append(item)

            # 4. Actions
            if "actions" in data:
                for act in data["actions"]:
                    self.add_action_row()
                    row = self.action_table.rowCount() - 1
                    
                    self.action_table.setItem(row, 0, QTableWidgetItem(act["time"]))
                    self.action_table.setItem(row, 2, QTableWidgetItem(act["state"]))
                    
                    combo = self.action_table.cellWidget(row, 1)
                    target_name = act["component"]
                    index = combo.findText(target_name)
                    if index >= 0:
                        combo.setCurrentIndex(index)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load: {str(e)}")
            import traceback
            traceback.print_exc()

    # --- SIMULATION & CONNECTIONS ---
    def run_simulation(self):
        sim_dialog = None
        if hasattr(self, "btn_run") and self.btn_run:
            self.btn_run.setEnabled(False)
        try:
            sim_dialog = self._show_sim_running_dialog()
            node_items = [i for i in self.scene.items() if isinstance(i, NodeItem)]
            conn_items = [i for i in self.scene.items() if isinstance(i, ConnectionLine)]
            
            if not node_items:
                return

            item_map = {} 
            conn_map = {}
            graph = {}
            
            # --- 1. INSTANTIATE REGULAR NODES (Skip Engine for now) ---
            regular_nodes = [i for i in node_items if i.node_type != "Engine"]
            engine_nodes = [i for i in node_items if i.node_type == "Engine"]
            
            for n_item in regular_nodes:
                cls = getattr(gfn, n_item.node_type)
                instance = cls(**n_item.params)
                item_map[n_item] = instance
                n_item.sim_instance = instance

            # --- 2. INSTANTIATE CONNECTIONS & BUILD GRAPH ---
            for c_item in conn_items:
                if c_item.conn_type == "Series":
                    sub_instances = []
                    for sub in c_item.params.get("connections", []):
                        sub_cls = getattr(gfn, sub["type"])
                        sub_instances.append(sub_cls(**sub["params"]))
                    
                    cls = getattr(gfn, "Series")
                    instance = cls(
                        connections=sub_instances, 
                        name=c_item.params.get("name", "series")
                    )
                else:
                    cls = getattr(gfn, c_item.conn_type)
                    instance = cls(**c_item.params)
                    
                c_item.sim_instance = instance
                conn_map[c_item] = instance

                # Map the graph ONLY if the ends are regular nodes. (Engines handle themselves)
                if c_item.start_item in item_map and c_item.end_item in item_map:
                    n1 = item_map[c_item.start_item]
                    n2 = item_map[c_item.end_item]
                    graph[instance] = (n1, n2)

            # --- 3. INSTANTIATE ENGINES (Injects Feed Lines dynamically) ---
            for eng_item in engine_nodes:
                # Find the connection lines drawn into the engine
                inbound_conns = [c for c in conn_items if c.end_item == eng_item]
                ox_conn_inst, fuel_conn_inst = None, None
                
                for c in inbound_conns:
                    source_item = c.start_item
                    # Try to guess which line is Ox and Fuel based on the source tank's fluid
                    fluid = source_item.params.get("fluid", "") or source_item.params.get("fluid_liq", "")
                    if fluid.lower() in ["oxygen", "lox", "n2o"]:
                        ox_conn_inst = conn_map[c]
                    else:
                        fuel_conn_inst = conn_map[c]
                
                # Instantiate Engine with injected connections
                instance = gfn.Engine(
                    fuel=eng_item.params.get("fuel", "n-Dodecane"),
                    oxidizer=eng_item.params.get("oxidizer", "Oxygen"),
                    ox_conn=ox_conn_inst,
                    fuel_conn=fuel_conn_inst,
                    eta_cstar=eng_item.params.get("eta_cstar", 0.92),
                    At=eng_item.params.get("At", 0.002),
                    Ae=eng_item.params.get("Ae", 0.01),
                    Pa=eng_item.params.get("Pa", 101325.0),
                    name=eng_item.params.get("name", "engine")
                )
                item_map[eng_item] = instance
                eng_item.sim_instance = instance
                
                # Route the inbound lines from the upstream tanks into the Engine node
                for c in inbound_conns:
                    n1 = item_map[c.start_item]
                    graph[conn_map[c]] = (n1, instance)

            # --- 4. PARSE ACTIONS (Builds the proper List of Tuples) ---
            actions = {}
            for row in range(self.action_table.rowCount()):
                t_str = self.action_table.item(row, 0).text()
                try: t_val = float(t_str)
                except: continue

                combo = self.action_table.cellWidget(row, 1)
                gui_item = combo.currentData()
                
                # Parse target state
                state_str = self.action_table.item(row, 2).text()
                try: state_val = float(state_str)
                except: state_val = 1.0 if state_str.lower() == "true" else 0.0
                
                sim_obj = None
                if isinstance(gui_item, tuple):
                    # It's a sub-component! Extract it from the parent Series
                    series_item, sub_idx = gui_item
                    if series_item in conn_map:
                        sim_obj = conn_map[series_item].connections[sub_idx]
                elif gui_item in conn_map:
                    # It's a top-level component
                    sim_obj = conn_map[gui_item]
                
                if sim_obj:
                    # THE FIX: Append to a LIST of tuples for the timestep!
                    if t_val not in actions:
                        actions[t_val] = []
                    actions[t_val].append((sim_obj, state_val))

            # Run
            net = gfn.Network(graph)
            t_max = float(self.inp_time.text())
            dt = float(self.inp_dt.text())
            net.sim(t_max, dt, actions=actions, verbose_steps=0)
            
            QMessageBox.information(self, "Success", "Simulation Complete!")

        except Exception as e:
            QMessageBox.critical(self, "Simulation Error", str(e))
            import traceback
            traceback.print_exc()
        finally:
            if sim_dialog:
                sim_dialog.close()
            if hasattr(self, "btn_run") and self.btn_run:
                self.btn_run.setEnabled(True)

    def create_connection(self, node1, node2):
        types = list(COMPONENT_DEFS["Connections"].keys())
        item, ok = QInputDialog.getItem(self, "Select Connection", "Type:", types, 0, False)
        if ok and item:
            defaults = COMPONENT_DEFS["Connections"][item].copy()
            conn_line = ConnectionLine(node1, node2, item, defaults)
            self.scene.addItem(conn_line)
            
            # Trigger Auto-Series Check
            self.check_and_form_series(node1)
            self.check_and_form_series(node2)

    def check_and_form_series(self, node):
        """Checks if a node is an intermediate pass-through and prompts to encapsulate into a Series."""
        if node.node_type != "Node":
            return
            
        attached = [item for item in self.scene.items() 
                    if isinstance(item, ConnectionLine) and (item.start_item == node or item.end_item == node)]
        
        if len(attached) == 2:
            c1, c2 = attached
            n_far_1 = c1.start_item if c1.end_item == node else c1.end_item
            n_far_2 = c2.start_item if c2.end_item == node else c2.end_item
            
            if n_far_1 == n_far_2: 
                return

            msg = QMessageBox.question(
                self, "Series Detected", 
                f"Node '{node.params.get('name', 'node')}' connects exactly two lines.\nWould you like to encapsulate them into a Series object?",
                QMessageBox.Yes | QMessageBox.No
            )
            
            if msg == QMessageBox.Yes:
                if c1.end_item == node:
                    first_c, second_c = c1, c2
                else:
                    first_c, second_c = c2, c1
                    
                series_conns = []
                for c in [first_c, second_c]:
                    if c.conn_type == "Series":
                        series_conns.extend(c.params.get("connections", []))
                    else:
                        series_conns.append({"type": c.conn_type, "params": c.params.copy()})
                        
                self.scene.removeItem(c1)
                self.scene.removeItem(c2)
                self.scene.removeItem(node)
                
                defaults = COMPONENT_DEFS["Connections"]["Series"].copy()
                defaults["connections"] = series_conns
                series_line = ConnectionLine(n_far_1, n_far_2, "Series", defaults)
                self.scene.addItem(series_line)
                
                if self.current_editing_item in (c1, c2, node):
                    self.current_editing_item = None
                    for i in reversed(range(self.prop_layout.count())): 
                        self.prop_layout.itemAt(i).widget().setParent(None)

    # --- PROPERTY EDITOR LOGIC ---
    def on_selection_changed(self):
        for i in reversed(range(self.prop_layout.count())): 
            self.prop_layout.itemAt(i).widget().setParent(None)
        
        items = self.scene.selectedItems()
        if not items:
            self.current_editing_item = None
            return

        item = items[0]
        self.current_editing_item = item
        
        lbl = f"<b>Type: {item.node_type}</b>" if isinstance(item, NodeItem) else f"<b>Type: {item.conn_type}</b>"
        self.prop_layout.addRow(QLabel(lbl))

        self.param_inputs = {}

        def render_params(params_dict, prefix=""):
            for key, val in params_dict.items():
                if key == "connections": continue 
                
                if key in ("fluid", "fluid_ullage"):
                    inp = QComboBox()
                    inp.setEditable(True)
                    inp.addItems(COMMON_FLUIDS)
                    inp.setCurrentText(str(val))
                    inp.currentTextChanged.connect(lambda text, k=key, d=params_dict: self.update_param(k, text, d))
                else:
                    inp = QLineEdit(str(val))
                    inp.textChanged.connect(lambda text, k=key, d=params_dict: self.update_param(k, text, d))
                
                label_text = _label_with_units(key)
                if prefix: label_text = f"{prefix} {label_text}"
                self.prop_layout.addRow(label_text, inp)
                
                input_key = f"{prefix}_{key}" if prefix else key
                self.param_inputs[input_key] = inp

        render_params(item.params)

        if isinstance(item, ConnectionLine) and item.conn_type == "Series":
            for idx, sub_conn in enumerate(item.params.get("connections", [])):
                lbl = QLabel(f"<b>--- Sub {idx+1}: {sub_conn['type']} ---</b>")
                lbl.setStyleSheet("color: #2980b9; margin-top: 10px;")
                self.prop_layout.addRow(lbl)
                render_params(sub_conn['params'], prefix=f"[{idx+1}]")

        if isinstance(item, NodeItem):
            self._prop_calc_item = item
            fluid_default = (
                item.params.get("fluid")
                or item.params.get("fluid_liq")
                or item.params.get("fluid_ullage")
                or "Water"
            )
            if self._prop_calc_widgets.get("fluid"):
                self._prop_calc_widgets["fluid"].setCurrentText(str(fluid_default))
        else:
            self._prop_calc_item = None

    def update_prop_calc_labels(self):
        w = self._prop_calc_widgets
        if not w:
            return
        w["lbl_p1_val"].setText(_label_with_units(w["p1"].currentText()))
        w["lbl_p2_val"].setText(_label_with_units(w["p2"].currentText()))
        w["lbl_out_val"].setText(_label_with_units(w["out"].currentText()))

    def compute_fluid_property(self):
        w = self._prop_calc_widgets
        if not w:
            return
        try:
            fluid = w["fluid"].currentText().strip()
            p1 = w["p1"].currentText().strip()
            p2 = w["p2"].currentText().strip()
            out = w["out"].currentText().strip()
            v1 = float(w["p1_val"].text().strip())
            v2 = float(w["p2_val"].text().strip())

            val = _try_props_refprop_then_coolprop(fluid, out, p1, v1, p2, v2)
            w["out_val"].setText(f"{val:.6g}")

            self._prop_calc_last.update({"p1": p1, "p2": p2, "out": out})

            if self._prop_calc_item and out in self._prop_calc_item.params:
                self._prop_calc_item.params[out] = val
                if out in self.param_inputs:
                    self.param_inputs[out].setText(str(val))
        except Exception as e:
            QMessageBox.warning(self, "Property Calc Error", str(e))

    def update_param(self, key, text, target_dict=None):
        if not self.current_editing_item:
            return
            
        if target_dict is None:
            target_dict = self.current_editing_item.params
            
        try: val = float(text)
        except ValueError: val = text
        
        target_dict[key] = val
        
        if key == "name" and isinstance(self.current_editing_item, NodeItem):
            self.current_editing_item.update()

    # --- UI KEY BINDS ---
    def keyPressEvent(self, event):
        if event.key() in (Qt.Key_Delete, Qt.Key_Backspace):
            self.delete_selected()
        else:
            super().keyPressEvent(event)

    def delete_selected(self):
        selected_items = self.scene.selectedItems()
        if not selected_items:
            return

        nodes_to_delete = set()
        conns_to_delete = set()

        for item in selected_items:
            if isinstance(item, NodeItem):
                nodes_to_delete.add(item)
            elif isinstance(item, ConnectionLine):
                conns_to_delete.add(item)

        all_conns = [i for i in self.scene.items() if isinstance(i, ConnectionLine)]
        
        for conn in all_conns:
            if conn in conns_to_delete: 
                continue 
            
            if conn.start_item in nodes_to_delete or conn.end_item in nodes_to_delete:
                conns_to_delete.add(conn)

        for conn in conns_to_delete:
            self.scene.removeItem(conn)
            
        for node in nodes_to_delete:
            self.scene.removeItem(node)
            
        if self.current_editing_item in nodes_to_delete or self.current_editing_item in conns_to_delete:
            self.current_editing_item = None
            for i in reversed(range(self.prop_layout.count())): 
                self.prop_layout.itemAt(i).widget().setParent(None)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainApp()
    window.show()
    sys.exit(app.exec_())
