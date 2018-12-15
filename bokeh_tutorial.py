"""Bokeh Visualization Template

This template is a general outline for turning your data into a 
visualization using Bokeh.
"""
# Data handling
#%%
import pandas as pd
import numpy as np

# Bokeh libraries
# main
from bokeh.io import output_file, output_notebook
from bokeh.plotting import figure, show
#secondary
from bokeh.models import ColumnDataSource
from bokeh.layouts import row, column, gridplot
from bokeh.models.widgets import Tabs, Panel

# reset output between sebsequent runs
# Import reset_output (only needed once) 
# from bokeh.plotting import reset_output
# Use reset_output() between subsequent show() calls, as needed
# reset_output()

# Prepare the data

# Determine where the visualization will be rendered (shows both if both are called)
#output_file('filename.html')  # Render to static HTML, and/or: 
#%%
output_notebook()  # Render inline in a Jupyter Notebook

# Set up the figure(s)
# fig = figure()  # Instantiate a figure() object
# example figure with custom settings:
# #%%
# fig = figure(background_fill_color='gray',
#              background_fill_alpha=0.5,
#              border_fill_color='blue',
#              border_fill_alpha=0.25,
#              plot_height=300,
#              plot_width=500,
#              h_symmetry=True,
#              x_axis_label='X Label',
#              x_axis_type='datetime',
#              x_axis_location='below',
#              x_range=('2018-01-01', '2018-06-30'),
#              y_axis_label='Y Label',
#              y_axis_type='linear',
#              y_axis_location='left',
#              y_range=(0, 100),
#              title='Example Figure',
#              title_location='above',
#              toolbar_location='below',
#              tools='save')

# change figure after instantiation:
# fig.grid.grid_line_color = None

# #%%
# # My x-y coordinate data
# x = [1, 2, 1]
# y = [1, 1, 2]

# # Output the visualization directly in the notebook
# output_file('first_glyphs.html', title='First Glyphs')

# # Create a figure with no toolbar and axis ranges of [0,3]
# fig = figure(title='My Coordinates',
#              plot_height=300, plot_width=300,
#              x_range=(0, 3), y_range=(0, 3),
#              toolbar_location=None)

# # Create first Glyphs
# # Draw the coordinates as circles
# fig.circle(x=x, y=y,
#            color='green', size=10, alpha=0.5)
#%%
# another example
# My word count data
day_num = np.linspace(1, 10, 10)
daily_words = [450, 628, 488, 210, 287, 791, 508, 639, 397, 943]
cumulative_words = np.cumsum(daily_words)

# Output the visualization directly in the notebook
output_notebook()

# Create a figure with a datetime type x-axis
fig = figure(title='My Tutorial Progress',
             plot_height=400, plot_width=700,
             x_axis_label='Day Number', y_axis_label='Words Written',
             x_minor_ticks=2, y_range=(0, 6000),
             toolbar_location=None)

# The daily words will be represented as vertical bars (columns)
fig.vbar(x=day_num, bottom=0, top=daily_words, 
         color='blue', width=0.75, 
         legend='Daily')

# The cumulative sum will be a trend line
fig.line(x=day_num, y=cumulative_words, 
         color='gray', line_width=1,
         legend='Cumulative')

# Put the legend in the upper left corner
fig.legend.location = 'top_left'


# Connect to and draw the data
# Organize the layout
# Preview and save 
show(fig)  # See what I made, and save if I like it