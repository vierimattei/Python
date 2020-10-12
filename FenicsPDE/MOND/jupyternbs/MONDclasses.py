#making a class for runtimes of each section to profile the program. it has a time and a name
#attribute, so when plotting we can directly use the name in e.g. a bar chart or pie chart
class BVP:
    
    #initialising class
    def __init__(self, weak_form, initial_guess, source, name):
        self.weak_form = weak_form
        self.initial_guess = initial_guess
        self.source = source
        self.name = name

#making a class for runtimes of each section to profile the program. it has a time and a name
#attribute, so when plotting we can directly use the name in e.g. a bar chart or pie chart
class run_time:
    
    #initialising class
    def __init__(self, time, name):
        self.time = time
        self.name = name
