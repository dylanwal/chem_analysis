


class MenuBarBase:
    def __init__(self, menubar):
        menu_file = menubar.addMenu("File")
        menu_processing = menubar.addMenu("Processing")
        menu_analysis = menubar.addMenu("Analysis")

        action_open = menu_file.addAction("Open", self.open)
        action_save = menu_file.addAction("Save")
        action_quit = menu_file.addAction("Quit", self.destroy)



        # Toolbars
        TBfile = window.addToolBar("File")
        TBfile.addAction(action_open)
        TBfile.addAction(action_save)

        open_icon = self.style().standardIcon(QtWidgets.QStyle.StandardPixmap.SP_DirOpenIcon)
        save_icon = self.style().standardIcon(QtWidgets.QStyle.StandardPixmap.SP_DialogSaveButton)

        action_open.setIcon(open_icon)
        action_save.setIcon(save_icon)

        TBfile.setAllowedAreas(QtCore.Qt.ToolBarArea.AllToolBarAreas)

        TBviewerNavigation = window.addToolBar("Viewer Navigation")

        action_previous = QtGui.QAction("<", self)
        action_next = QtGui.QAction(">", self)

        TBviewerNavigation.addAction(action_previous)
        TBviewerNavigation.addAction(action_next)
        TBviewerNavigation.setAllowedAreas(
            QtCore.Qt.ToolBarArea.TopToolBarArea | QtCore.Qt.ToolBarArea.BottomToolBarArea
        )