import matplotlib.legend_handler
from matplotlib.container import ErrorbarContainer
import matplotlib as mat
# mat.use("TkAgg")
import matplotlib.pyplot as plt
# %matplotlib inline
plt.switch_backend("TkAgg")
# plt.switch_backend("Qt5Agg")


mat.rcParams.update({'font.size': 16})
mat.rcParams["font.family"] = "Times New Roman"
size = 12

# see this link to see why this works: https://stackoverflow.com/questions/19470104/python-matplotlib-errorbar-legend-picking


class re_order_errorbarHandler(matplotlib.legend_handler.HandlerErrorbar):
    """
    Sub-class the standard error-bar handler 
    """

    def create_artists(self, *args, **kwargs):
        #  call the parent class function
        a_list = matplotlib.legend_handler.HandlerErrorbar.create_artists(
            self, *args, **kwargs)
        # re-order the artist list, only the first artist is added to the
        # legend artist list, this is the one that corresponds to the markers
        a_list = a_list[-1:] + a_list[:-1]
        return a_list


# def errorbar(x, d, err=None, fmt='o-', capthick=1, capsize=4, label=None, size=4, shift=False, ax=plt.gca(), **kwargs):
#     if err is None:
#         err = d*0.0
#     ax.errorbar(x, d, yerr=err, fmt='o-',
#                 capthick=1, capsize=4, color=color, **kwargs)
#     return ax


color = ['k', 'r', 'b', 'g', 'm', 'c', 'navy',
         'y', 'cyan', 'darkgreen', 'violet', 'lime', 'purple']
color = color*40


class InteractiveLegend(object):
    """click the label to hide a errorbar curve"""

    def __init__(self, ax):
        handles, labels = ax.get_legend_handles_labels()
        # each element of handles are the returned obj from plt.errobar(...)
        leg = ax.get_legend()
        leg.legendHandles.append(
            {ErrorbarContainer: re_order_errorbarHandler(numpoints=2)})
        self.fig = leg.axes.figure
        self.lined = dict()
        for legline, origline in zip(leg.legendHandles, handles):
            legline.set_picker(5)  # 5 pts tolerance
            self.lined[legline] = origline
        # print self.lined

        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        # self.fig.canvas.draw()

    def onpick(self, event):
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        legline = event.artist
        # print legline
        origline = self.lined[legline]
        for a in origline.get_children():
            vis = not a.get_visible()
            a.set_visible(vis)
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        self.fig.canvas.draw()
