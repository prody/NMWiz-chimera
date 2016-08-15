# -----------------------------------------------------------------------------
#
from chimera.extension import EMO, manager

# -----------------------------------------------------------------------------
#
class NMWiz_EMO(EMO):

    def name(self):
        return 'Normal Mode Wizard'
    def description(self):
        return 'Visualization of Normal Mode'
    def categories(self):
        return ['Structure Analysis']
    def icon(self):
        return None
    def activate(self):
        self.module('gui').show_nmwiz_dialog()
        return None

# -----------------------------------------------------------------------------
#
manager.registerExtension(NMWiz_EMO(__file__))
