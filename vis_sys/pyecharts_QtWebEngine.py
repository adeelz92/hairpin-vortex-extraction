"""
In this version of checkpoint, we extract a subtree from
the mega tree and make that visible only. We specify the
minimum size of a vortex and then apply the filtration.
We move the PCP plot down in this version.
Integrated new clustering results in this version.
In this version, we make the plotting adaptive to the number of attributes in the vortex profile.
In this version, we added functionality of selecting multiple vortices.
In this version, we added the functionality to view only the selected vortices.
In this version, we changed the structure of tree json and added value item in every node.
In this version, we changed the self.vortexProfiles from list to dict.
In this version, we sorted the tree node with the currently selected option properly.
"""

import sys
import numpy as np
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtWebEngineWidgets import QWebEngineView, QWebEnginePage

from pyecharts.globals import CurrentConfig
from pyecharts.charts import Tree, Parallel, Scatter
from pyecharts import options as opts

# noinspection PyUnresolvedReferences
import vtk
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonCore import vtkCollectionIterator
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera
from vtkmodules.vtkFiltersModeling import vtkOutlineFilter
from vtkmodules.vtkFiltersCore import (
    vtkConnectivityFilter,
    vtkPolyDataConnectivityFilter,
    vtkCleanPolyData,
    vtkThreshold,
)
from vtkmodules.vtkRenderingAnnotation import vtkCubeAxesActor
from vtkmodules.vtkIOLegacy import vtkDataSetReader, vtkCompositeDataReader
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet, vtkPolyData
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkActorCollection,
    vtkPolyDataMapper,
    vtkDataSetMapper,
    vtkCellPicker,
    vtkProperty,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkRenderer
)

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5 import Qt
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import os
import codecs
import json
from operator import itemgetter


# CurrentConfig.ONLINE_HOST = "https://fastly.jsdelivr.net/npm/echarts@5.4.0/dist/"


class _LoggedPage(QWebEnginePage):

    def __init__(self, parent=None):
        self.parent = parent
        super(_LoggedPage, self).__init__()

    def javaScriptConsoleMessage(self, level: 'QWebEnginePage.JavaScriptConsoleMessageLevel', message: str,
                                 lineNumber: int, sourceID: str) -> None:
        print(message)
        if "_" in message:
            # regionId, level, parent
            regionId, level, parent = message.strip().split("_")[0:3]
            self.parent.update_dataset(int(regionId), int(level), int(parent), flag="tree")
        elif "," in message:
            # regionId, level, parent
            regionId, level, parent = message.strip().split(",")[0:3]
            self.parent.update_dataset(int(regionId), int(level), int(parent), flag="pcp")
        else:
            idx = int(message.strip())
            clusters = self.parent.getVortexProfile(idx)
            self.parent.update_dataset(clusters=clusters, flag="cluster")
            # print(regionId, level, parent)


class MyQVTKRenderWindowInteractor(QVTKRenderWindowInteractor):
    def __init__(self, ren, pick_actor_func, remove_actor_func):
        super(MyQVTKRenderWindowInteractor, self).__init__()
        self.ren = ren
        self.pick_actor_func = pick_actor_func
        self.remove_actor_func = remove_actor_func

    def mouseDoubleClickEvent(self, event: QEvent):
        if event.button() == 1:
            clickPos = self.GetEventPosition()
            picker = vtkCellPicker()
            picker.Pick(clickPos[0], clickPos[1], 0, self.ren)
            cellId = picker.GetCellId()
            if picker.GetCellId() != -1:
                self.pick_actor_func(cellId)
        elif event.button() == 2:
            clickPos = self.GetEventPosition()
            picker = vtkCellPicker()
            picker.Pick(clickPos[0], clickPos[1], 0, self.ren)
            cellId = picker.GetCellId()
            if picker.GetCellId() != -1:
                self.remove_actor_func(cellId)


class PlotCanvas(FigureCanvas):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, fig)
        self.axes = self.figure.add_subplot(111)
        self.axes.set_title('Vortex Profile')
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        # self.plot()

    def plot(self, data):
        normalized_data = data / np.linalg.norm(data)

        self.axes.clear()
        self.axes.set_title('Vortex Profile')
        self.axes.boxplot(normalized_data, sym="")
        self.draw()


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        Qt.QMainWindow.__init__(self, parent)

        ''' Step 1: Create main layout '''
        self.setWindowTitle('Regions Tree')
        # self.resize(1000, 1000)
        self.frame = Qt.QFrame()  # Create a main window frame to add ui widgets
        self.mainLayout = Qt.QHBoxLayout()  # Set layout - Lines up widgets horizontally
        self.frame.setLayout(self.mainLayout)
        self.setCentralWidget(self.frame)

        ''' Step 2: Create left layout for dataset and tree visualization '''
        self.leftWidget = Qt.QWidget()
        self.left_panel_layout = Qt.QVBoxLayout()
        self.leftWidget.setLayout(self.left_panel_layout)
        self.mainLayout.addWidget(self.leftWidget, stretch=3)

        ''' Step 3: Add a vtk widget to the left widget '''
        # As we use QVBoxLayout, the vtk widget will be automatically moved to the top
        self.ren1 = vtkRenderer()
        self.ren2 = vtkRenderer()
        self.vtkWidget = MyQVTKRenderWindowInteractor(self.ren1, self.pick_regions, self.remove_regions)
        ''' To do: Configure self.vtkWidget for self.ren2'''
        self.left_panel_layout.addWidget(self.vtkWidget, stretch=1)

        ''' Step 4: Add a graph view widget to the central widget '''
        self.webPlotView = QWebEngineView()
        self.plotPage = _LoggedPage(self)
        self.webPlotView.setPage(self.plotPage)
        self.left_panel_layout.addWidget(self.webPlotView, stretch=1)

        # Initialize the vtk variables for the visualization
        self.init_vtk_widget()

        ''' Step 5: Create right layout for controls and graphs visualization '''
        self.rightWidget = Qt.QWidget()
        self.right_panel_layout = Qt.QVBoxLayout()
        self.rightWidget.setLayout(self.right_panel_layout)
        self.mainLayout.addWidget(self.rightWidget, stretch=2)

        ''' Step 6: Add controls to the interface '''
        self.add_controls()

        ''' Step 7: Add a web view widget to the left widget '''
        self.webTreeView = QWebEngineView()
        self.treePage = _LoggedPage(self)
        self.webTreeView.setPage(self.treePage)
        self.right_panel_layout.addWidget(self.webTreeView, stretch=2)

        ''' Step 8: Add a graph view widget to the central widget '''
        # self.graph = PlotCanvas()
        self.webScatterView = QWebEngineView()
        self.scatterPage = _LoggedPage(self)
        self.webScatterView.setPage(self.scatterPage)
        self.right_panel_layout.addWidget(self.webScatterView, stretch=2)

    def init_vtk_widget(self):
        vtk.vtkObject.GlobalWarningDisplayOff()  # Disable vtkOutputWindow

        # Create the graphics structure. The renderer renders into the render
        # window. The render window interactor captures mouse events and will
        # perform appropriate camera or actor manipulation depending on the
        # nature of the events.
        self.colors = vtkNamedColors()
        self.dataset = None
        self.selection_mode = "single"
        self.view_mode = "all"
        self.render_mode = "grid"

        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren1)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()

        # The following set the interactor for 2D image style (i.e., no rotation)
        style = vtkInteractorStyleTrackballCamera()
        # style = MouseInteractorHighLightActor(self.ren1, self.pick_regions)
        self.iren.SetInteractorStyle(style)
        self.ren1.SetBackground(self.colors.GetColor3d("White"))  # you can change the background color here
        self.ren2.SetBackground(self.colors.GetColor3d("White"))

        # Start the vtk screen
        self.ren1.ResetCamera()
        self.ren2.ResetCamera()
        self.iren.Initialize()
        self.iren.Start()

        # self.init_dataset()

    def add_controls(self):
        ''' Add a sample group box '''
        groupBox = Qt.QGroupBox("Options")  # Use a group box to group controls
        groupBox_layout = Qt.QVBoxLayout()  # lines up the controls vertically
        groupBox_layout.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop)
        groupBox.setLayout(groupBox_layout)
        self.right_panel_layout.addWidget(groupBox, stretch=1)

        ''' Add a textfield ( QLineEdit) to show the file path and the browser button '''
        openFileBox_layout = Qt.QHBoxLayout()
        label = Qt.QLabel("Choose a file (e.g., vtk):")
        openFileBox_layout.addWidget(label, stretch=1)
        self.qt_file_name = Qt.QLineEdit()
        openFileBox_layout.addWidget(self.qt_file_name, stretch=8)
        qt_browser_button = Qt.QPushButton('Browse')
        qt_browser_button.clicked.connect(self.on_file_browser_clicked)
        qt_browser_button.show()
        openFileBox_layout.addWidget(qt_browser_button, stretch=1)
        ''' Add the Open button '''
        qt_open_button = Qt.QPushButton('Open')
        qt_open_button.clicked.connect(self.open_vtk_file)
        qt_open_button.show()
        openFileBox_layout.addWidget(qt_open_button, stretch=1)
        file_widget = Qt.QWidget()
        file_widget.setLayout(openFileBox_layout)
        groupBox_layout.addWidget(file_widget)

        ''' Create a new widget for sort menu'''
        sortBox_layout = Qt.QHBoxLayout()
        sortBoxlabel = Qt.QLabel("Sort By: ")
        sortBoxlabel.setMaximumWidth(int(sortBoxlabel.width() / 15))
        sortBox_layout.addWidget(sortBoxlabel, stretch=1)

        '''
        36 criteria are being used for vortex profile construction which are below
        0. lambda_2				1. lambda_ci			2. Q
        3. delta				4. divergence			5. oyf
        6. size(# of voxels)	7. mag(vorticity)		8. enstrophy
        9. mag(Velocity)		10. mag(Acceleration) 	11. norm(Jacobian)
        12. Avg(Velocity_x)		13. Avg(Velocity_y)		14. Avg(Velocity_z)
        15. Avg(Accelera_x)		16. Avg(Accelera_y)		17. Avg(Accelera_z)
        18. Avg(vorticity_x)	19. Avg(vorticity_y)	20. Avg(vorticity_z)
        21. Avg(Jacobian_11)	22. Avg(Jacobian_12)	23. Avg(Jacobian_13)
        24. Avg(Jacobian_21)	25. Avg(Jacobian_22)	26. Avg(Jacobian_23)
        27. Avg(Jacobian_31)	28. Avg(Jacobian_32)	29. Avg(Jacobian_33)
        30. Min(lambda_2)		31. Max(lambda_2)		32. Min(Q)
        33. Max(Q)				34. Min(oyf)			35. Max(oyf)
        '''

        self.itemsList = []

        self.sortMenu = Qt.QComboBox()
        self.sortMenu.addItems(self.itemsList)
        sortBox_layout.addWidget(self.sortMenu, stretch=1)

        nNodeslabel = Qt.QLabel("Number of nodes: ")
        # nNodeslabel.setMargin(int(sortBoxlabel.width() / 10))
        sortBox_layout.addWidget(nNodeslabel, stretch=1)

        self.nNodes_scale = Qt.QDoubleSpinBox()
        # set the initial values of some parameters
        self.nNodes_scale.setValue(10)
        self.nNodes_scale.setRange(1, 100)
        self.nNodes_scale.setSingleStep(1)
        sortBox_layout.addWidget(self.nNodes_scale, stretch=1)

        sizelabel = Qt.QLabel("Specify min size: ")
        # sizelabel.setMargin(int(sortBoxlabel.width() / 10))
        sortBox_layout.addWidget(sizelabel, stretch=1)

        self.size_scale = Qt.QDoubleSpinBox()
        # set the initial values of some parameters
        self.size_scale.setRange(50, 30000)
        self.size_scale.setValue(500)
        self.size_scale.setSingleStep(50)
        sortBox_layout.addWidget(self.size_scale, stretch=1)

        ''' Add the Open button '''
        qt_sort_button = Qt.QPushButton('Apply')
        qt_sort_button.setFixedWidth(qt_open_button.width())
        qt_sort_button.clicked.connect(self.sort_nodes)
        qt_sort_button.show()
        sortBox_layout.addWidget(qt_sort_button, stretch=1)
        sortWidget = Qt.QWidget()
        sortWidget.setLayout(sortBox_layout)
        groupBox_layout.addWidget(sortWidget)

        ''' Add View and Selection Mode layouts'''
        viewAndSelect_layout = Qt.QHBoxLayout()
        # Add the view mode layout
        viewMode_layout = Qt.QGridLayout()
        viewMode_label = Qt.QLabel("View: ")
        viewAndSelect_layout.addWidget(viewMode_label, stretch=1)
        # Add the first radio button
        radioButton = Qt.QRadioButton("Selected")
        radioButton.setChecked(False)
        radioButton.mode = "selected"
        radioButton.toggled.connect(self.on_view_mode_changed)
        viewMode_layout.addWidget(radioButton, 0, 0)
        # Add the second radio button
        radioButton = Qt.QRadioButton("All")
        radioButton.setChecked(True)
        radioButton.mode = "all"
        radioButton.toggled.connect(self.on_view_mode_changed)
        viewMode_layout.addWidget(radioButton, 0, 1)
        viewMode_widget = Qt.QFrame()
        viewMode_widget.setFrameStyle(Qt.QFrame.Panel | Qt.QFrame.Raised)
        viewMode_widget.setLayout(viewMode_layout)
        viewAndSelect_layout.addWidget(viewMode_widget, stretch=2)

        # Add the render mode layout
        renderMode_layout = Qt.QGridLayout()
        renderMode_label = Qt.QLabel("Render: ")
        viewAndSelect_layout.addWidget(renderMode_label, stretch=1)
        # Add the first radio button
        radioButton = Qt.QRadioButton("Grid")
        radioButton.setChecked(True)
        radioButton.mode = "grid"
        radioButton.toggled.connect(self.on_render_mode_changed)
        renderMode_layout.addWidget(radioButton, 0, 0)
        # Add the second radio button
        radioButton = Qt.QRadioButton("Streamlines")
        radioButton.setChecked(False)
        radioButton.mode = "streamlines"
        radioButton.toggled.connect(self.on_render_mode_changed)
        renderMode_layout.addWidget(radioButton, 0, 1)
        # Add the third radio button
        radioButton = Qt.QRadioButton("Isosurface")
        radioButton.setChecked(False)
        radioButton.mode = "isosurface"
        radioButton.toggled.connect(self.on_render_mode_changed)
        renderMode_layout.addWidget(radioButton, 0, 2)
        self.renderMode_widget = Qt.QFrame()
        self.renderMode_widget.setFrameStyle(Qt.QFrame.Panel | Qt.QFrame.Raised)
        self.renderMode_widget.setLayout(renderMode_layout)
        self.renderMode_widget.setDisabled(True)
        viewAndSelect_layout.addWidget(self.renderMode_widget, stretch=2)

        # Add the selection mode layout
        selectMode_layout = Qt.QGridLayout()
        selectMode_label = Qt.QLabel("Select: ")
        viewAndSelect_layout.addWidget(selectMode_label, stretch=1)
        # Add the first radio button
        radioButton = Qt.QRadioButton("Single")
        radioButton.setChecked(True)
        radioButton.mode = "single"
        radioButton.toggled.connect(self.on_selection_mode_changed)
        selectMode_layout.addWidget(radioButton, 0, 0)
        # Add the second radio button
        radioButton = Qt.QRadioButton("Multiple")
        radioButton.setChecked(False)
        radioButton.mode = "multiple"
        radioButton.toggled.connect(self.on_selection_mode_changed)
        selectMode_layout.addWidget(radioButton, 0, 1)
        selectMode_widget = Qt.QFrame()
        selectMode_widget.setFrameStyle(Qt.QFrame.Panel | Qt.QFrame.Raised)
        selectMode_widget.setLayout(selectMode_layout)
        viewAndSelect_layout.addWidget(selectMode_widget, stretch=2)
        # Add the clear button
        clear_selection_button = Qt.QPushButton('Clear')
        clear_selection_button.setFixedWidth(qt_open_button.width())
        clear_selection_button.clicked.connect(self.clear_selection)
        clear_selection_button.show()
        viewAndSelect_layout.addWidget(clear_selection_button, stretch=1)
        viewAndSelect_widget = Qt.QWidget()
        viewAndSelect_widget.setLayout(viewAndSelect_layout)
        groupBox_layout.addWidget(viewAndSelect_widget)

        # Show the window
        self.showMaximized()

    def on_file_browser_clicked(self):
        dlg = Qt.QFileDialog()
        dlg.setFileMode(Qt.QFileDialog.AnyFile)
        dlg.setNameFilter("loadable files (*.vtk *.mhd)")

        if dlg.exec_():
            filenames = dlg.selectedFiles()
            self.qt_file_name.setText(filenames[0])

    def getItemList(self, data):
        name = data['name']
        if name == "root":
            children = data['children']
            if children != None or children != []:
                for child in children:
                    # Recursion
                    return self.getItemList(child)
        else:
            itemlist = []
            data_list = str(name).strip().split("_")
            for value in data_list:
                if "=" in value:
                    criteria_name = value.strip().split("=")[0]
                    if criteria_name != "Isovalue":
                        itemlist.append(criteria_name)
            return itemlist

    def open_vtk_file(self):
        '''Read and verify the vtk input file '''
        input_file_name = self.qt_file_name.text()
        json_file_name = input_file_name[0:-4] + ".json"
        clusters_file_name = input_file_name[0:-4] + "_clusters.txt"
        self.reader = vtkCompositeDataReader()
        self.reader.SetFileName(input_file_name)
        self.reader.Update()
        self.multi_block_dataset = vtkMultiBlockDataSet.SafeDownCast(self.reader.GetOutput())

        # Remove all actors if already exist
        ren1Actors = self.ren1.GetActors()
        for ren1Actor in ren1Actors:
            self.ren1.RemoveActor(ren1Actor)
        ren2Actors = self.ren2.GetActors()
        for ren2Actor in ren2Actors:
            self.ren2.RemoveActor(ren2Actor)

        # Render the dataset
        self.init_dataset()

        with codecs.open(json_file_name, "r", encoding="utf-8") as f:
            j = json.load(f)

        self.childList = []
        # Adds options to the menu
        self.reviseTreeJson(j)
        self.treeData = [j]
        print(self.childList)
        self.sortMenu.addItems(self.itemsList)

        # Find the index of Lambda2
        self.size_index = -1
        for index, item in enumerate(self.itemsList):
            if item == "Size":
                self.size_index = index
                break

        # Build vortex criteria
        self.vortexProfiles = dict()
        self.buildVortexProfiles(self.treeData[0])
        clusters = np.loadtxt(clusters_file_name)
        self.cluster_data = [list(x) for x in clusters]

        # Sort and load the tree
        self.sort_nodes()

        # Load the scatter plot
        clusters = [x[-1] for x in self.cluster_data]
        elements, counts = np.unique(clusters, return_index=True)
        nClusters = len(counts)
        clus_data = [arr[-3:] for arr in self.cluster_data]
        self.update_and_load_scatter(data=clus_data, n_clusters=nClusters)

    def on_selection_mode_changed(self):
        radioButton = self.sender()
        if radioButton.isChecked() == True:
            self.selection_mode = radioButton.mode
            print(self.selection_mode)

    def on_view_mode_changed(self):
        radioButton = self.sender()
        if radioButton.isChecked() == True:
            self.view_mode = radioButton.mode
            if self.view_mode == "selected":
                self.vtkWidget.GetRenderWindow().RemoveRenderer(self.ren1)
                self.vtkWidget.GetRenderWindow().AddRenderer(self.ren2)
                self.renderMode_widget.setDisabled(False)
                # self.ren2.ResetCamera()
            elif self.view_mode == "all":
                self.vtkWidget.GetRenderWindow().RemoveRenderer(self.ren2)
                self.vtkWidget.GetRenderWindow().AddRenderer(self.ren1)
                self.renderMode_widget.setDisabled(True)
            self.vtkWidget.GetRenderWindow().Render()
            # (self.view_mode)


    def get_unique_actors(self):
        all_actors = []
        if hasattr(self, 'previos_actor'):
            if self.previos_actor is not None:
                all_actors.append(self.previos_actor)
        if hasattr(self, 'selected_actors'):
            if len(self.selected_actors) > 0:
                for actor in self.selected_actors:
                    if actor not in all_actors:
                        all_actors.append(actor)
        return all_actors


    def on_render_mode_changed(self):
        radioButton = self.sender()
        if radioButton.isChecked() == True:
            self.render_mode = radioButton.mode
            if hasattr(self, 'previos_actor'):
                if self.previos_actor is not None:
                    data = self.previos_actor.GetMapper().GetInput()
                    print(data)

    def clear_selection(self):
        for _, actor in self.selected_actors.items():
            self.ren1.RemoveActor(actor)
            self.ren2.RemoveActor(actor)
        self.selected_actors.clear()
        self.ren1.RemoveActor(self.previos_actor)
        self.ren2.RemoveActor(self.previos_actor)
        self.previos_actor = None
        self.vtkWidget.GetRenderWindow().Render()
        self.plotPCPdata(regionIds=())

    def update_and_load_chart(self, data):
        textStyle = dict()
        textStyle.update(fontWeight='bold')

        width = self.webPlotView.width()
        height = self.webPlotView.height()
        parallelChart = Parallel(init_opts=opts.InitOpts(animation_opts=opts.AnimationOpts(animation=False),
                                                         width="{}px".format(width), height="{}px".format(height)))
        parallelChart.options.update(parallel=opts.ParallelOpts(pos_top="15%", pos_right="5%"))

        numdims = len(data[0]["value"])

        axis_labels = opts.LabelOpts(font_size=8)
        axis_labels.update(width=25, overflow="truncate")
        textStyle = {"fontWeight": "bold", "fontSize": 12}
        # textStyle.update(nameTextStyle={})
        pcp = (
            parallelChart
                .add_schema(
                schema=[dict({'dim':i, 'name':str(self.itemsList[i-3]), 'type':'value',
                              'nameLocation': 'end' if i%2 ==0 else 'start',
                              'axisLabel': axis_labels, 'nameTextStyle':textStyle})
                        for i in range(3, numdims)])
                .add("", data, linestyle_opts=opts.LineStyleOpts(opacity=0.5))
                .set_global_opts(title_opts=opts.TitleOpts(title="PCP View"),
                                 legend_opts=opts.LegendOpts(is_show=False))
        )

        """
        pcp = (
            parallelChart
                .add_schema(
                [
                    dict({'dim': 3, 'name': "\u03BB_2", 'type':'value', 'nameLocation': 'end', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 4, 'name': "\u03BB_ci", 'type':'value', 'nameLocation': 'start', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 5, 'name': "Q", 'type':'value', 'nameLocation': 'end', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 6, 'name': "\u0394", 'type':'value', 'nameLocation': 'start', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 7, 'name': "Div", 'type':'value', 'nameLocation': 'end', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 8, 'name': "\u03C9yf", 'type':'value', 'nameLocation': 'start', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 9, 'name': "Size", 'type':'value', 'nameLocation': 'end', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 10, 'name': "\u03C9", 'type':'value', 'nameLocation': 'start', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 11, 'name': "\u03BE", 'type':'value', 'nameLocation': 'end', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 12, 'name': "V", 'type':'value', 'nameLocation': 'start', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 13, 'name': "a", 'type':'value', 'nameLocation': 'end', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 14, 'name': "J", 'type':'value', 'nameLocation': 'start', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 15, 'name': "Vx", 'type':'value', 'nameLocation': 'end', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 16, 'name': "Vy", 'type':'value', 'nameLocation': 'start', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 17, 'name': "Vz", 'type':'value', 'nameLocation': 'end', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 18, 'name': "Ax", 'type':'value', 'nameLocation': 'start', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 19, 'name': "Ay", 'type':'value', 'nameLocation': 'end', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 20, 'name': "Az", 'type':'value', 'nameLocation': 'start', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 21, 'name': "\u03C9x", 'type':'value', 'nameLocation': 'end', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 22, 'name': "\u03C9y", 'type':'value', 'nameLocation': 'start', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 23, 'name': "\u03C9z", 'type':'value', 'nameLocation': 'end', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 24, 'name': "J11", 'type':'value', 'nameLocation': 'start', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 25, 'name': "J12", 'type':'value', 'nameLocation': 'end', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 26, 'name': "J13", 'type':'value', 'nameLocation': 'start', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 27, 'name': "J21", 'type':'value', 'nameLocation': 'end', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 28, 'name': "J22", 'type':'value', 'nameLocation': 'start', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 29, 'name': "J23", 'type':'value', 'nameLocation': 'end', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 30, 'name': "J31", 'type':'value', 'nameLocation': 'start', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 31, 'name': "J32", 'type':'value', 'nameLocation': 'end', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                    dict({'dim': 32, 'name': "J33", 'type':'value', 'nameLocation': 'start', 'axisLabel': a, 'axisTick': b, 'nameTextStyle':textStyle}),
                ],
            )
                .add("", data, linestyle_opts=opts.LineStyleOpts(opacity=0.5))
                .set_global_opts(title_opts=opts.TitleOpts(title="PCP View"),
                                 legend_opts=opts.LegendOpts(is_show=False))
        )
        """

        # fmt: off
        pcp.add_js_funcs("chart_" + pcp.chart_id + ".on('click', function(params) { console.log(params.value); });")
        # fmt: on
        # Render the line chart and save as html file
        html_file = (os.getcwd() + "\\pcp.html")
        pcp.render(path=html_file)

        # Insert into the QtWebEngine window
        # self.webTreeView.load(QUrl.fromLocalFile(html_file))

        self.plotPage.load(QUrl.fromLocalFile(html_file))
        self.webPlotView.update()

    def update_and_load_tree(self, data):

        # print(self.webTreeView.size())
        width = self.webTreeView.width()
        height = self.webTreeView.height()

        tree = (
            Tree(init_opts=opts.InitOpts(width="{}px".format(width), height="{}px".format(height)))
                .add("", data,
                     # pos_top="-10%", pos_bottom="-10%",
                     symbol_size=15, is_roam=True,
                     is_expand_and_collapse=False,
                     label_opts=opts.LabelOpts(is_show=False),
                     leaves_label_opts=opts.LabelOpts(is_show=False),
                     itemstyle_opts=opts.ItemStyleOpts(color='Red'))
                .set_global_opts(title_opts=opts.TitleOpts(title="Tree View"))
        )
        # tree.options['series'][0]['height'] = "90%"
        # fmt: off
        tree.add_js_funcs("chart_" + tree.chart_id + ".on('click', function(params) { console.log(params.value); });")
        # fmt: on
        # Render the line chart and save as html file
        html_file = (os.getcwd() + "\\tree.html")
        tree.render(path=html_file)

        # Insert into the QtWebEngine window
        # self.webTreeView.load(QUrl.fromLocalFile(html_file))
        self.treePage.load(QUrl.fromLocalFile(html_file))
        self.webTreeView.update()

    def update_and_load_scatter(self, data, n_clusters):

        if data == None:
            data = list(self.cluster_data)
        # print(data)
        '''
        data = [
            [10.0, 8.04, 1],
            [8.0, 6.95, 1],
            [13.0, 7.58, 1],
            [9.0, 8.81, 2],
            [11.0, 8.33, 2],
            [14.0, 9.96, 2],
            [6.0, 7.24, 3],
            [4.0, 4.26, 3],
            [12.0, 10.84, 3],
            [7.0, 4.82, 4],
            [5.0, 5.68, 4],
        ]
        '''
        # data.sort(key=lambda x: x[0])
        # print(self.webTreeView.size())
        width = int(self.webScatterView.width())
        height = int(self.webScatterView.height())

        scatter = Scatter(init_opts=opts.InitOpts(width="{}px".format(width), height="{}px".format(height)))
        scatter.add_dataset(source=data)
        scatter.add_yaxis(series_name="", y_axis=data, symbol_size=15,
                          label_opts=opts.LabelOpts(is_show=False),
                          itemstyle_opts=opts.ItemStyleOpts(border_color='#555')
                          )
        scatter.options['xAxis'][0]['scale'] = True
        scatter.options['yAxis'][0]['scale'] = True
        pieces = []
        for i in range(-1, n_clusters - 1):
            dd = dict()
            dd.update(value=i)
            dd.update(label='cluster ' + str(i))
            dd.update(color=scatter.colors[i])
            pieces.append(dd)

        scatter.set_global_opts(title_opts=opts.TitleOpts(title="Scatter Plot"),
                                visualmap_opts=opts.VisualMapOpts(type_="color",
                                                                  is_piecewise=True,
                                                                  # min_=0,
                                                                  # max_=n_clusters - 2,
                                                                  pos_top='middle',
                                                                  pos_left='5',
                                                                  # split_number=n_clusters,
                                                                  pieces=pieces,
                                                                  dimension=2
                                                                  ))
        # Move the grid to the right
        scatter.options.update(grid=opts.GridOpts(pos_left=150))

        # fmt: off
        scatter.add_js_funcs(
            "chart_" + scatter.chart_id + ".on('click', function(params) { console.log(params.dataIndex); });")
        # fmt: on
        # Render the line chart and save as html file
        html_file = (os.getcwd() + "\\scatter.html")
        scatter.render(path=html_file)

        # Insert into the QtWebEngine window
        # self.webTreeView.load(QUrl.fromLocalFile(html_file))
        self.scatterPage.load(QUrl.fromLocalFile(html_file))
        self.webScatterView.update()

    def getVortexProfile(self, idx):
        class_label = self.cluster_data[idx][-1]
        print("Class label: ", class_label)
        final_profiles = [cluster[0:3] for cluster in self.cluster_data if cluster[-1] == class_label]
        total_elements = len(self.cluster_data[0][3:-3])
        min_list = [np.inf for i in range(total_elements)]
        max_list = [-np.inf for i in range(total_elements)]
        for cluster in self.cluster_data:
            if cluster[-1] == class_label:
                for j, ele in enumerate(cluster[3:-3]):
                    if min_list[j] > ele:
                        min_list[j] = ele
                    if max_list[j] < ele:
                        max_list[j] = ele
        print(min_list)
        print(max_list)
        print(10 * "#")
        # a, b, c = self.cluster_data[idx][0:3]
        # a, b, c = self.vortexProfiles[idx][-4:-1]
        return final_profiles

    def init_dataset(self):

        self.selected_actors = dict()
        self.previos_actor = None
        self.PCPdata = None
        self.dataset =  vtkPolyData.SafeDownCast(self.multi_block_dataset.GetBlock(1))

        datasetMapper = vtkDataSetMapper()
        datasetMapper.SetInputData(self.dataset)
        datasetMapper.ScalarVisibilityOff()
        # datasetMapper.SetScalarModeToUseCellData()
        # datasetMapper.SetScalarRange(scalar_range[0], scalar_range[1])

        datasetActor = vtkActor()
        datasetActor.SetMapper(datasetMapper)
        datasetActor.GetProperty().SetColor(self.colors.GetColor3d("Grey"))
        datasetActor.GetProperty().SetOpacity(0.15)

        outline = vtkOutlineFilter()
        outline.SetInputData(self.dataset)
        outline.Update()
        outlineMapper = vtkPolyDataMapper()
        outlineMapper.SetInputConnection(outline.GetOutputPort())
        outlineActor = vtkActor()
        outlineActor.SetMapper(outlineMapper)
        outlineActor.GetProperty().SetColor(self.colors.GetColor3d("Black"))
        outlineActor.GetProperty().SetLineWidth(2.0)

        cubeAxesActor = vtkCubeAxesActor()
        cubeAxesActor.SetUseTextActor3D(0)
        cubeAxesActor.SetBounds(outline.GetOutput().GetBounds())
        cubeAxesActor.SetLabelOffset(10)
        cubeAxesActor.SetTitleOffset(10)
        cubeAxesActor.SetLabelScaling(False, 0, 0, 0)
        cubeAxesActor.GetTitleTextProperty(0).SetColor(self.colors.GetColor3d("Black"))
        cubeAxesActor.GetLabelTextProperty(0).SetColor(self.colors.GetColor3d("Black"))
        cubeAxesActor.GetTitleTextProperty(1).SetColor(self.colors.GetColor3d("Black"))
        cubeAxesActor.GetLabelTextProperty(1).SetColor(self.colors.GetColor3d("Black"))
        cubeAxesActor.GetTitleTextProperty(2).SetColor(self.colors.GetColor3d("Black"))
        cubeAxesActor.GetLabelTextProperty(2).SetColor(self.colors.GetColor3d("Black"))
        cubeAxesActor.XAxisMinorTickVisibilityOff()
        cubeAxesActor.YAxisMinorTickVisibilityOff()
        cubeAxesActor.ZAxisMinorTickVisibilityOff()
        # cubeAxesActor.SetScreenSize(10.0);
        cubeAxesActor.SetFlyMode(vtkCubeAxesActor.VTK_FLY_OUTER_EDGES)

        cubeAxesActor.SetCamera(self.ren1.GetActiveCamera())
        self.ren1.AddActor(datasetActor)
        self.ren1.AddActor(outlineActor)
        self.ren1.AddActor(cubeAxesActor)
        self.ren1.GetActiveCamera().SetPosition(0, -1, 0)
        self.ren1.GetActiveCamera().SetViewUp(0, 0, 1)
        self.ren1.ResetCamera()
        self.ren1.GetActiveCamera().Azimuth(-30)
        self.ren1.GetActiveCamera().Elevation(15)
        self.ren1.GetActiveCamera().Zoom(1.5)

        # Re-render the screen
        if self.view_mode == "selected":
            self.ren2.ResetCamera()
        # self.vtkWidget.GetRenderWindow().Render()

    def getActor(self, level=1, regionId=0, color="Blue"):
        arr_name = "RegionIds_{}".format(level)

        """
        if level == 0:
            arr_name = "RegionIds"
        """

        self.temp_dataset = vtkPolyData.SafeDownCast(self.multi_block_dataset.GetBlock(level))
        threshold = vtkThreshold()
        threshold.SetInputData(self.temp_dataset)
        threshold.SetInputArrayToProcess(0, 0, 0, 1, arr_name)
        threshold.SetLowerThreshold(regionId)
        threshold.SetUpperThreshold(regionId)
        threshold.Update()

        data = threshold.GetOutput()

        mapper = vtkDataSetMapper()
        mapper.SetInputData(data)
        mapper.ScalarVisibilityOff()

        actor = vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(self.colors.GetColor3d(color))
        actor.GetProperty().SetOpacity(0.9)
        return actor

    def removeActor(self, regionId):
        if regionId in self.selected_actors.keys():
            actor = self.selected_actors.pop(regionId)
            self.ren1.RemoveActor(actor)
            self.ren2.RemoveActor(actor)

    def pick_regions(self, cellId):
        if self.previos_actor != None:
            # self.previos_actor.GetProperty().SetColor(self.colors.GetColor3d("Grey"))
            self.ren1.RemoveActor(self.previos_actor)
            self.ren2.RemoveActor(self.previos_actor)
            self.previos_actor = None

        if self.selection_mode == "single":
            for regionId, actor in self.selected_actors.items():
                self.ren1.RemoveActor(actor)
                self.ren2.RemoveActor(actor)
            self.selected_actors.clear()

        arr_name = "RegionIds_1"
        colors_array = self.dataset.GetCellData().GetArray(arr_name)
        regionId = colors_array.GetTuple1(cellId)
        if regionId == 0:
            return
        actor = self.getActor(level=1, regionId=regionId)
        self.selected_actors[regionId] = actor
        self.ren1.AddActor(actor)
        self.ren2.AddActor(actor)

        # Re-render the screen
        if self.view_mode == "selected":
            self.ren2.ResetCamera()
        self.vtkWidget.GetRenderWindow().Render()
        regionIds = self.selected_actors.keys()
        self.plotPCPdata(tuple(regionIds))

    def remove_regions(self, cellId):
        arr_name = "RegionIds"
        colors_array = self.dataset.GetCellData().GetArray(arr_name)
        regionId = colors_array.GetTuple1(cellId)
        if regionId == 0:
            return
        self.removeActor(regionId)

        # Re-render the screen
        if self.view_mode == "selected":
            self.ren2.ResetCamera()
        self.vtkWidget.GetRenderWindow().Render()
        regionIds = self.selected_actors.keys()
        self.plotPCPdata(tuple(regionIds))

    def update_dataset(self, regionId=0, level=0, parent=0, clusters=None, flag=None):

        if (clusters != None):
            for _, actor in self.selected_actors.items():
                self.ren1.RemoveActor(actor)
                self.ren2.RemoveActor(actor)
            self.selected_actors.clear()

        def getScalarAvgValue(array):
            endId = array.GetNumberOfTuples()
            TotalValue = 0
            for tupleId in range(endId):
                thisValue = array.GetTuple1(tupleId)
                TotalValue += thisValue
            TotalValue /= array.GetNumberOfTuples()
            return TotalValue

        # print(regionId, level, parent)
        if (level > 0):
            if self.previos_actor != None:
                if len(self.selected_actors) == 0:
                    self.ren1.RemoveActor(self.previos_actor)
                    self.ren2.RemoveActor(self.previos_actor)
                    # self.previos_actor.GetProperty().SetColor(self.colors.GetColor3d("Grey"))
                else:
                    selected = False
                    for _, a in self.selected_actors.items():
                        if a == self.previos_actor:
                            selected = True
                            break
                    if selected == False:
                        self.ren1.RemoveActor(self.previos_actor)
                        self.ren2.RemoveActor(self.previos_actor)
                    else:
                        self.previos_actor.GetProperty().SetColor(self.colors.GetColor3d("Blue"))

            if regionId in self.selected_actors.keys() and flag == "pcp":
                current_actor = self.selected_actors[regionId]
                actors = self.ren1.GetActors()
                for actor in actors:
                    if actor == current_actor:
                        actor.GetProperty().SetColor(self.colors.GetColor3d("Red"))
                        self.previos_actor = actor
                self.updatePCPdata(regionId)
            else:
                actor = self.getActor(int(level), int(regionId), color="Red")
                # self.selected_actors[regionId] = actor
                self.ren1.AddActor(actor)
                self.ren2.AddActor(actor)
                selected = False
                for _, a in self.selected_actors.items():
                    if a == self.previos_actor:
                        selected = True
                        break
                if selected == False:
                    self.ren1.RemoveActor(self.previos_actor)
                    self.ren2.RemoveActor(self.previos_actor)

                self.previos_actor = actor
                self.updatePCPdata(regionId)

        elif clusters != None:
            self.webScatterView.update()
            cluster_reg_ids = []
            topIds = []
            for cluster in clusters:
                regionId, level, parent = cluster
                id = str(int(regionId)) + "_" + str(int(level)) + "_" + str(int(parent))
                topIds.append(id)
                actor = self.getActor(int(level), int(regionId))
                self.selected_actors[regionId] = actor
                self.ren1.AddActor(actor)
                self.ren2.AddActor(actor)
                cluster_reg_ids.append(int(regionId))
            '''To do'''
            ''' Plot PCP for all clusters in the data'''
            # Get the sorted subtree
            newSubTree = None
            tempTreeData = self.treeData.copy()
            for id in topIds:
                newSubTree = self.getSubTree(tempTreeData[0], id, newSubTree)
            self.updateTreeJson(newSubTree, topIds)
            self.update_and_load_tree([newSubTree])
            self.plotPCPdata(cluster_reg_ids)


        # Re-render the screen
        self.vtkWidget.GetRenderWindow().Render()
        self.webTreeView.update()

    def buildVortexProfiles(self, data):

        for element in data:
            if element == "name":
                name = data[element]
                if name == "root":
                    continue
                # Get all values
                data_list = data['value']
                ids = data_list[0:4]
                data_values = data_list[4:]
                profile = []
                for value in data_values:
                    profile.append(float(value))
                for id in ids:
                    profile.append(id)
                # Append the profile to the final dataset
                self.vortexProfiles[ids[0]] = profile
            elif element == "children":
                children = data[element]
                if children is not None:
                    for child in data[element]:
                        # Recursion
                        self.buildVortexProfiles(child)

    def plotPCPdata(self, regionIds):
        length = 0
        for profile in self.vortexProfiles:
            length = len(self.vortexProfiles[profile]) - 1
            break

        if len(regionIds) == 0:
            forPCP = [{"value": [-1 for i in range(length)],
                       "lineStyle": opts.LineStyleOpts(color="rgba(0, 0, 255, 1)", width=1, opacity=0.5)}]
        else:
            forPCP = [{"value": self.vortexProfiles[regionId][-4:-1] + self.vortexProfiles[regionId][0:-4],
                       "lineStyle": opts.LineStyleOpts(color="rgba(0, 0, 255, 1)", width=1, opacity=0.5)}
                      for regionId in regionIds]

        self.PCPdata = forPCP
        self.update_and_load_chart(self.PCPdata)

    def addPCPdata(self, regionId):
        if len(self.selected_actors.keys()) == 0:
            pcp = []
        else:
            pcp = self.PCPdata.copy()
        for key, currprofile in self.vortexProfiles.items():
            if currprofile[-4] == regionId:
                pcpDict = dict()
                pcpDict.update(value=currprofile[-4:-1] + currprofile[0:-4])
                pcpDict.update(lineStyle=opts.LineStyleOpts(color="rgba(255, 0, 0, 1)", width=2, opacity=1))
                pcp.append(pcpDict)
                break
        return pcp

    def updatePCPdata(self, regionId):
        found = False
        for i, pcpDict in enumerate(self.PCPdata):
            pcpDict.update(lineStyle=opts.LineStyleOpts(color="rgba(0, 0, 255, 1)", width=1, opacity=0.5))
            # values = pcpDict["value"]
            # print(values, regionId)
            if pcpDict["value"][0] == regionId:
                pcpDict.update(lineStyle=opts.LineStyleOpts(color="rgba(255, 0, 0, 1)", width=2, opacity=1))
                found = True
            self.PCPdata[i] = pcpDict

        if found == False:
            pcp = self.addPCPdata(regionId=regionId)
            self.update_and_load_chart(pcp)
            return
        # Replot the PCP chart
        self.update_and_load_chart(self.PCPdata)

    # Find top nodes specified by the number of Nodes menu
    def findTopNodes(self, n, criteriaIdx):

        # Make a temporary copy of vortex profiles
        # tempProfiles = self.vortexProfiles.copy()
        # First find all those profiles which has size > min Size

        tempVortexProfiles = [profile for key, profile in self.vortexProfiles.items() if
                              profile[self.size_index] > self.size_scale.value()]

        '''
        for profile in self.vortexProfiles:
            if profile[6] > self.size_scale.value():
                tempVortexProfiles.append(profile)
        '''

        topIds = []

        # Sort the vortex profiles in descending order to find top n nodes
        sortedList = sorted(tempVortexProfiles, key=itemgetter(criteriaIdx), reverse=True)
        # tempVortexProfiles.sort(key = lambda tempVortexProfile : tempVortexProfiles[criteriaIdx])

        forPCP = []
        count = 0
        for currprofile in sortedList:
            if (int(currprofile[-4]) not in self.childList):
                continue

            pcpDict = dict()
            pcpDict.update(value=currprofile[-4:-1] + currprofile[0:-4])
            pcpDict.update(lineStyle=opts.LineStyleOpts(color="rgba(0, 0, 255, 1)", width=1, opacity=0.5))
            forPCP.append(pcpDict)

            id = str(int(currprofile[-4])) + "_" + str(int(currprofile[-3])) + "_" + str(int(currprofile[-2]))
            topIds.append(id)
            count += 1
            if (count == n):
                break
        # print("****************")
        return topIds, forPCP

    def reviseTreeJson(self, data: dict):
        self.itemsList = self.getItemList(data)

        name = data['name']
        if name == "root":
            data.update(value=None)
            children = data['children']
            if children != None or children != []:
                for child in children:
                    # Recursion
                    self.reviseTreeJson(child)
            return
        else:
            data_list = str(name).strip().split("_")
            profile = []
            for value in data_list:
                if "=" in value:
                    key_val = value.strip().split("=")
                    val = float(key_val[1])
                    profile.append(val)
                else:
                    profile.append(int(value))

            data.update(name=profile[0])
            data.update(value=profile)
            children = data['children']
            if children is not None and children != []:
                for child in children:
                    # Recursion
                    self.reviseTreeJson(child)
            else:
                self.childList.append(profile[0])

    def updateTreeJson(self, data, topIds):
        name = data['name']

        if name == "root":
            children = data['children']
            if children is not None:
                for child in children:
                    # Recursion
                    self.updateTreeJson(child, topIds)
            return

        # Get all values
        data_list = data['value']
        ids = data_list[0:3]

        id = str(ids[0]) + "_" + str(ids[1]) + "_" + str(ids[2])
        if id not in topIds:
            # print(data)
            data.update(itemStyle=opts.series_options.ItemStyleOpts(opacity=1))
            data.update(collapsed=False)
        else:
            data.update(itemStyle=opts.series_options.ItemStyleOpts(opacity=1))
            data.update(collapsed=False)

        children = data['children']
        if children is not None:
            for child in children:
                # Recursion
                self.updateTreeJson(child, topIds)

    def updateTree(self, newTree: dict, chosenList: list = None, finalIdx=None, currentIdx=0):
        if finalIdx == None:
            for i in range(len(chosenList) - 1, -1, -1):
                if chosenList[i]["chosen"] == True:
                    finalIdx = i + 1
                    break

        if finalIdx == None or currentIdx == finalIdx:
            return
        else:
            name = chosenList[currentIdx]["name"]
            values = chosenList[currentIdx]["value"]
            # print(len(newTree["children"]))
            # Check if child already exists
            found = False
            for i in range(len(newTree["children"])):
                if newTree["children"][i]["name"] == name:
                    found = True
                    break

            currentIdx += 1
            if found == False:
                child = dict()
                child.update(name=name)
                child.update(value=values)
                child.update(children=[])
                # Add a new child
                newTree["children"].append(child)
                self.updateTree(newTree["children"][-1], chosenList, finalIdx, currentIdx)
            else:
                self.updateTree(newTree["children"][i], chosenList, finalIdx, currentIdx)

    def getSubTree(self, data: dict, topId, newTree=None, chosenList=None):
        if newTree == None:
            newTree = dict()
            newTree.update(name="root")
            newTree.update(value=None)
            newTree.update(children=[])

        # print(newTree)
        if chosenList == None:
            chosenList = []

        name = data['name']
        if name == "root":
            children = data['children']
            if children is not None:
                for child in children:
                    # Recursion
                    self.getSubTree(child, topId, newTree, chosenList)
            return newTree

        # Get all values
        data_list = data['value']
        selected = dict()
        selected.update(name=name)
        selected.update(value=data_list)
        ids = data_list[0:3]
        id = str(ids[0]) + "_" + str(ids[1]) + "_" + str(ids[2])
        # regionId, level, parent = data_list[0:3]
        # id = str(int(regionId)) + "_" + str(int(level)) + "_" + str(int(parent))
        found = False
        if id != topId:
            # data.update(itemStyle=opts.series_options.ItemStyleOpts(opacity=0))
            # data.update(collapsed=True)
            selected.update(chosen=False)
        else:
            # data.update(itemStyle=opts.series_options.ItemStyleOpts(opacity=1))
            # data.update(collapsed=False)
            found = True
            selected.update(chosen=True)
        chosenList.append(selected)

        children = data['children']
        if children is not None:
            for child in children:
                # Recursion
                self.getSubTree(child, topId, newTree, chosenList)
        if found:
            self.updateTree(newTree, chosenList)
            return

        chosenList.remove(selected)

    def sort_nodes(self):

        self.clear_selection()

        '''
        # Clear all actors
        if self.previos_actor != None:
                self.previos_actor.GetProperty().SetColor(self.colors.GetColor3d("Grey"))
        for regionId, actor in self.selected_actors.items():
            self.ren1.RemoveActor(actor)
        self.selected_actors.clear()
        '''

        # Make a temporary copy of the tree data
        tempTreeData = self.treeData.copy()
        selectedIdx = self.sortMenu.currentIndex()
        nNodes = int(self.nNodes_scale.value())

        topIds, self.PCPdata = self.findTopNodes(nNodes, selectedIdx)
        # print(topIds)
        # Get the sorted subtree
        newSubTree = None
        for id in topIds:
            newSubTree = self.getSubTree(tempTreeData[0], id, newSubTree)
        self.updateTreeJson(newSubTree, topIds)
        self.update_and_load_tree([newSubTree])
        self.update_and_load_chart(self.PCPdata)
        '''
        if self.view_mode == "selected":
            self.ren2.ResetCamera()
        '''
        self.vtkWidget.GetRenderWindow().Render()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
