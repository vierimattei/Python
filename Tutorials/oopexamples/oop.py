# -*- coding: utf-8 -*-

class Employee:
    
    num_of_emps = 0
    raise_amt = 1.04

    def __init__(self, first, last):
        self.first = first
        self.last = last
        # self.pay = pay
        
        #increasing num_of_emps every time we create an instance of the class
        #(an object)
        # Employee.num_of_emps += 1
      
    @property
    def email(self):
        return '{}.{}@email.com'.format(self.first, self.last)  
      
    @property
    def fullname(self):
        return '{} {}'.format(self.first, self.last)
    
    @fullname.setter
    def fullname(self,name):
        first, last = name.split(' ')
        self.first = first
        self.last = last
        
    @fullname.deleter
    def fullname(self):
        print('Delete Name')
        self.first = None
        self.last = None
    
    # def apply_raise(self):
    #     self.pay = int(self.pay * self.raise_amt)
        
    # #methods (functions inside a class) that have __ before and after are
    # #called 'dunder' (cause double under score) methods. They are 'special' 
    # #because they are connected to a special keyword or symbol in the 
    # #language. e.g. repr to change how we print our object to give more
    # #meaningful information to a developer
    # def __repr__(self):
    #     return "Employee('{}', '{}', '{}')".format(self.first, self.last, self.pay)
    
    # #str is also for printing but to give better information to a user
    # def __str__(self):
    #     return "{} - {}".format(self.fullname(), self.email)

    # #defining how to add two employee objects together (arbitrary, here pay)
    # #same dunder __add__ is used when we e.g. add two strings
    # def __add__(self, other):
    #     return self.pay + other.pay
            
    # #here saying that the length of the employee object is the # characters in
    # #the full name
    # def __len__ (self):
    #     return len(self.fullname())

emp_1 = Employee('John','Smith')

emp_1.first = 'Jim'

emp_1.fullname = 'Corey Schafer'

print(emp_1.first)
print(emp_1.email)
print(emp_1.fullname)

del emp_1.fullname

# print(emp_1+emp_2)

# print(issubclass(Manager, Developer))

# print(emp_1)

# print(repr(emp_1))
# print(str(emp_2))

# print(int.__add__(1,2))
# print(str.__add__('a','b'))

# @classmethod
# def set_raise_amt(cls, amount):
#     cls.raise_amt = amount
    
# @classmethod
# def from_string(cls, emp_str):
#     first, last, pay = emp_str.split('-')
#     return cls(first, last, pay)

# @staticmethod
# def is_workday(day):
#     if day.weekday() == 5 or day.weekday() == 6:
#         return False
#     else:
#         return True

#Subclass inheriting from Employee
# class Developer(Employee):
#     raise_amt = 1.1
    
#     def __init__(self, first, last, pay, prog_lang):
#         super().__init__(first, last, pay)
#         self.prog_lang = prog_lang
        
# #Other subclass
# class Manager(Employee):
    
#     def __init__(self, first, last, pay, employees = None):
#         super().__init__(first, last, pay)
#         if employees == None:
#             self.employees = []
#         else:
#             self.employees = employees
    
#     def add_emp(self, emp):
#         if emp not in self.employees:
#             self.employees.append(emp)
            
#     def remove_emp(self, emp):
#         if emp in self.employees:
#             self.employees.remove(emp)
            
#     def print_emps(self):
#         for emp in self.employees:
#             print('-->', emp.fullname())

# print(mgr_1.email)

# mgr_1.add_emp(dev_2)
# mgr_1.remove_emp(dev_1)

# mgr_1.print_emps()

# print(dev_1.email)
# print(dev_1.prog_lang)

# print(dev_1.pay)
# dev_1.apply_raise()
# print(dev_1.pay)


# import datetime
# my_date = datetime.date(2016, 7, 11)

# print(Employee.is_workday(my_date))


# emp_str_1 = 'John-Doe-70000'
# emp_str_2 = 'Steve-Smith-30000'
# emp_str_3 = 'Jane-Doe-90000'

# new_emp_1 = Employee.from_string(emp_str_1)

# Employee.set_raise_amt(1.05)

# print(Employee.raise_amt)
# print(dev_1.email)
# print(dev_2.email)

# print(emp_1.__dict__)

