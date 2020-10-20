#making a class for runtimes of each section to profile the program. it has a time and a name
#attribute, so when plotting we can directly use the name in e.g. a bar chart or pie chart
class BVP:
    
    #inisitalisation
    def __init__(self, weak_form, initial_guess, source, name):
        self.weak_form = weak_form
        self.initial_guess = initial_guess
        self.source = source
        self.name = name

#making a class for runtimes of each section to profile the program. it has a time and a name
#attribute, so when plotting we can directly use the name in e.g. a bar chart or pie chart
class run_time:
    
    #inisitalisation
    def __init__(self, time, name):
        self.time = time
        self.name = name

#Making a class for each quantity we want to plot. We want it to have a numpy array for its values
#and a string we can use for the subplot title or label in a plot. Also color for plotting, and line
#width.
class plot_quantity:
    
    #inisitalisation
    def __init__(self, quantity, title, color, width):
        self.quantity = quantity
        self.title = title
        self.color = color
        self.width = width