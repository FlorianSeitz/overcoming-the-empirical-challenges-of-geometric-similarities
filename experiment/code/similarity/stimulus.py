#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on Jan 26, 2016

@author: ralbrecht
'''


from expyriment import control, stimuli, misc
from cmath import sqrt
from numpy.ma.core import cos, sin
import matplotlib
import pygame

darkgreen = pygame.Color(matplotlib.colors.rgb2hex([0,.5,0]))
darkred = pygame.Color(matplotlib.colors.rgb2hex([.5,0,0]))
grid_width = 35
grid_height = 35
num_cueVals = 5
instructions_ok = "OK"
instructions_feedback_wrong  = "Wrong"
instructions_feedback_own_answer = "Hallo"
instructions_feedback_correct_answer = "Correct Answer"
instructions_next = "Next"
instructions_feedback_correct = "Correct"

class AdjustableStimulus:
    def __init__(self, side, offset_y = 75, offset_x=100, pos=1, forced_choice=False):
        self.docu = []
        self.my_control = 0
        self.canvas_stim, self.response_rect = self.get_canvas(pos, offset_x, offset_y, forced_choice, side)
        self.u = grid_width
        self.h = abs(sqrt(3)/2*self.u)
        self.cells1 = []
        self.cells2 = []
        self.cells3 = []
        self.cells4 = []
        self.shapes1 = []
        self.shapes2 = []
        self.shapes3 = []
        self.shapes4 = []
        self.grid = []
        self.size1 = 0
        self.size2 = 0
        self.size3 = 0
        self.size4 = 0
        self.size_corr1 = 0
        self.size_corr2 = 0
        self.size_corr3 = 0
        self.size_corr4 = 0
        self.button_ok = stimuli.Rectangle(size=(100,50), position=(0,-300), colour=misc.constants.C_GREY)
        self.ok_text = stimuli.TextLine(text=instructions_ok, position=self.button_ok.position, text_colour=misc.constants.C_BLACK)
    
    def get_canvas(self, pos, offset_x, offset_y, forced_choice, side):
        if side == "right":
            canvas_stim = stimuli.Canvas((control.defaults.window_size[0]/2, 428.5), (400, 75), colour=misc.constants.C_WHITE)
        if side == "left":
            canvas_stim = stimuli.Canvas((control.defaults.window_size[0]/2, 428.5), (-400, 75), colour=misc.constants.C_WHITE)
        if side == "middle":
            canvas_stim = stimuli.Canvas((control.defaults.window_size[0]/2, 428.5), (0, -50), colour=misc.constants.C_WHITE)
        response_rect = 0

        return (canvas_stim, response_rect)
        
    
        
    ### Create Shapes
    def create_shape_triangle(self, pos, rotate, color, height=grid_height, width=grid_width):
        shape = stimuli.Shape(position=(0,0), colour=color)
        shape.add_vertex((width/2, 1.732*height/2))
        shape.add_vertex((width/2, -1.732*height/2))
        shape.rotate(rotate)
        shape.reposition(pos)
#        shape.preload()
        return shape
    
    def create_shape_hex(self, pos, rotate, color, height=grid_height, width=grid_width):
        r=(height/2)
        shape = stimuli.Shape(position=(r*cos(0),r*sin(0)), colour=color)
        shape.add_vertex((r*cos((3.14/180)*72),r*sin((3.14/180)*72)))
        shape.add_vertex((r*cos((3.14/180)*144),r*sin((3.14/180)*144)))
        shape.add_vertex((r*cos((3.14/180)*216),r*sin((3.14/180)*216)))
        shape.add_vertex((r*cos((3.14/180)*288),r*sin((3.14/180)*288)))
        shape.rotate(rotate)
        shape.reposition(pos)
#        shape.preload()
        return shape
    
    def create_shape_cross(self, pos, rotate, color, height=grid_height, width=grid_width):
        shape = stimuli.FixCross((width + 1.5, height -2.5), pos, (grid_width+grid_height)/4, color)
        shape.rotate(rotate)
        shape.reposition(pos)
#        shape.preload()
        return shape
    
    def create_shape_rect(self, height, pos, rotate, color):
        shape = stimuli.Rectangle((height-5, height-5), position=(0,0), colour=color)
        shape.rotate(rotate)
        shape.reposition((pos[0], pos[1]))
#        shape.preload()
        return shape
    
    def create_shape_circ(self, height, pos, rotate, color):
        shape = stimuli.Circle(height/2, position=(0,0), colour=color)
        shape.reposition(pos)
#        shape.preload()
        return shape
    
    def create_cell(self, pos, rotate, col):
        r = stimuli.Rectangle((grid_width, grid_height), misc.constants.C_WHITE,  1, position=(0,0))
        r.rotate(rotate)
        r.reposition(new_position=pos)
#        r.preload()
        self.grid.append(r)
        r_in = stimuli.Rectangle((grid_width, grid_height), col,  0, position=(0,0))
        r_in.rotate(rotate)
        r_in.reposition(new_position=pos)
#        r_in.preload()
        return r_in
    
    
    def create_cues(self, deg, cells, offset =(0,0), col = (0, 0, 0, 255)):
        i = 0
        alpha = (3.14/180)*(deg)
        px = offset[0]+(grid_width*(cos(-alpha)))
        py = offset[1]+(grid_width*(sin(-alpha)))
        while (i < num_cueVals):
            px = px+(grid_width*(cos(-alpha)))
            py = py-(grid_width*(sin(-alpha)))
            r = self.create_cell((px,py), deg, col)
            cells.append(r)
            i = i + 1
                   
    
    def create_grid(self, pal_set, pos_set):
        # ADJUST
        # For rotation, better don't change!
        #Test
        pos1 = (157, 23, 296.3, 243.3) 
        # Vertical and horizontal position of arms and leg grid   
        pos2 = ((grid_width, grid_height),(-grid_width, grid_height), (-5,-grid_height-20),(7,-grid_height-20)) 
        # Position of shapes within grid. You should not need to change it
        pos3 = (157-90,270+23,296.3+90,243+90)
        
        pal = [pygame.Color("#c7813f"),pygame.Color("#3fa44a"), pygame.Color("#469bb4"), pygame.Color("#e94de3")]
        col_1 = pal[pal_set[0]]
        col_2 = pal[pal_set[1]]
        col_3 = pal[pal_set[2]]
        col_4 = pal[pal_set[3]]
        self.create_cues(pos1[pos_set[0]], self.cells1, pos2[pos_set[0]], pygame.Color(255 - col_1[0], 255 - col_1[1], 255 - col_1[2], col_1[3]))
        self.create_cues(pos1[pos_set[1]], self.cells2, pos2[pos_set[1]], pygame.Color(255 - col_2[0], 255 - col_2[1], 255 - col_2[2], col_2[3]))
        self.create_cues(pos1[pos_set[2]], self.cells3, pos2[pos_set[2]], pygame.Color(255 - col_3[0], 255 - col_3[1], 255 - col_3[2], col_3[3]))
        self.create_cues(pos1[pos_set[3]], self.cells4, pos2[pos_set[3]], pygame.Color(255 - col_4[0], 255 - col_4[1], 255 - col_4[2], col_4[3]))
        
        pos_docu = ("left-arm", "right-arm", "right-leg", "left-leg")
        pal_docu = ("orange", "green", "blue", "violet")
            
        self.docu=[[pos_docu[pos_set[0]], pal_docu[pal_set[0]], self.size1+1, "rectangle"],
        [pos_docu[pos_set[1]], pal_docu[pal_set[1]], self.size2+1, "circle"],
        [pos_docu[pos_set[2]], pal_docu[pal_set[2]], self.size3+1, "triangle"],
        [pos_docu[pos_set[3]], pal_docu[pal_set[3]], self.size4+1, "cross"]]
        
        for c in self.cells1:
            r = self.create_shape_rect(grid_height, c.position , pos3[pos_set[0]], col_1)
            self.shapes1.append(r)    
        for c in self.cells2:
            r = self.create_shape_rect(grid_height, c.position , pos3[pos_set[1]], col_2)
            self.shapes2.append(r)
        for c in self.cells3:
            r = self.create_shape_rect(grid_height, c.position , pos3[pos_set[2]], col_3)
            self.shapes3.append(r)
        for c in self.cells4:
            r = self.create_shape_cross(c.position , pos3[pos_set[3]]+90, col_4)
            self.shapes4.append(r)
    
    def create_grid_3(self, pal_set, pos_set):
        # ADJUST
        # For rotation, better don't change!
        #Test
        pos1 = (150, 30, 270)
        # Vertical and horizontal position of arms and leg grid   
        pos2 = ((18.5, 24.2),(-16.7, 24.2), (0,-grid_height-13.5))
        #pos2 = ((grid_width, grid_height),(-grid_width, grid_height), (0,-grid_height-15))
        # Position of shapes within grid. You should not need to change it
        pos3 = (150-90,270+30,0)
        
        pal = [pygame.Color("#c7813f"),pygame.Color("#3fa44a"), pygame.Color("#469bb4")]
        col_1 = pal[pal_set[0]]
        col_2 = pal[pal_set[1]]
        col_3 = pal[pal_set[2]]
        self.create_cues(pos1[pos_set[0]], self.cells1, pos2[pos_set[0]], pygame.Color(255 - col_1[0], 255 - col_1[1], 255 - col_1[2], col_1[3]-100))
        self.create_cues(pos1[pos_set[1]], self.cells2, pos2[pos_set[1]], pygame.Color(255 - col_2[0], 255 - col_2[1], 255 - col_2[2], col_2[3]-100))
        self.create_cues(pos1[pos_set[2]], self.cells3, pos2[pos_set[2]], pygame.Color(255 - col_3[0], 255 - col_3[1], 255 - col_3[2], col_3[3]-100))
        
        pos_docu = ("left-arm", "right-arm", "leg")
        pal_docu = ("orange", "green", "blue")
            
        self.docu=[[pos_docu[pos_set[0]], pal_docu[pal_set[0]], self.size1+1, "rectangle"],
        [pos_docu[pos_set[1]], pal_docu[pal_set[1]], self.size2+1, "circle"],
        [pos_docu[pos_set[2]], pal_docu[pal_set[2]], self.size3+1, "triangle"]]
        
        for c in self.cells1:
            r = self.create_shape_circ(grid_height, c.position, pos3[pos_set[0]], col_1)
            self.shapes1.append(r)    
        for c in self.cells2:
            r = self.create_shape_cross(c.position , pos3[pos_set[1]], col_2)
            self.shapes2.append(r)
        for c in self.cells3:
            r = self.create_shape_triangle(c.position , pos3[pos_set[2]], col_3, height=grid_height + 1.5, width=grid_width)
            self.shapes3.append(r)
    
    def show_grid(self, size1, size2, size3, bool_ok = False) :
        # Draw upper body
        triangleTop = self.create_shape_triangle((0,grid_height/2), 0, misc.constants.C_GREY, height=grid_height*3, width= grid_width*3)
#        triangleTop.preload()
        triangleTop.plot(self.canvas_stim)
        # Draw lower body (ADJUST)
        #triangleBottom = self.create_shape_triangle((0, -grid_height-5), 180, misc.constants.C_GREY, height= grid_height/2, width=grid_width*2.5)
        #triangleBottom.preload()
        #triangleBottom.plot(self.canvas_stim)
        # Draw head
        #        head = stimuli.Rectangle((grid_height*2, grid_height*2.2), position=(0,2.75*grid_height), colour=misc.constants.C_GREY)
        #        head.preload()
        #        head.plot(self.canvas_stim)
        #        eye1 = self.create_shape_circ(10, (-grid_width/4,3*grid_height), 0, misc.constants.C_BLACK)
        #        eye1.preload()
        #        eye1.plot(self.canvas_stim)
        #        eye2 = self.create_shape_circ(10, (grid_width/4,3*grid_height), 0, misc.constants.C_BLACK)
        #        eye2.preload()
        #        eye2.plot(self.canvas_stim)
        #        nose = self.create_shape_triangle((0, 2.65*grid_height), 180, misc.constants.C_BLACK, 10, 8)
        #        nose.preload()
        #        nose.plot(self.canvas_stim)
        #        mouth = stimuli.Ellipse((10,4), misc.constants.C_BLACK, position = (0, 2.3*grid_height))
        #        mouth.preload()
        #        mouth.plot(self.canvas_stim)

        if bool_ok:
            self.button_ok.plot(self.canvas_stim)
            self.ok_text.plot(self.canvas_stim)
        
        for r in self.cells1:
            r.plot(self.canvas_stim)
        for r in self.cells2:
            r.plot(self.canvas_stim)
        for r in self.cells3:
            r.plot(self.canvas_stim)
        for r in self.grid:
            r.plot(self.canvas_stim)
        
        for s in range(size1+1):
            self.shapes1[s].plot(self.canvas_stim)
        for s in range(size2+1):
            self.shapes2[s].plot(self.canvas_stim)
        for s in range(size3+1):
            self.shapes3[s].plot(self.canvas_stim)
        
        
    # TODO: Only relevant for instructions?
    def give_feedback(self, c_s1, c_s2, c_s3, c_s4, text, color, instruction):
        
        l_button = stimuli.Rectangle(size=(300,50), position=(0,150), colour=misc.constants.C_GREY)
        l_text = stimuli.TextLine(text=text, position=l_button.position, text_colour=misc.constants.C_BLACK)
        
        f_button = stimuli.Rectangle(size=(300,50), position=(0,200), colour=misc.constants.C_GREY)
        f_text = stimuli.TextLine(text=instructions_feedback_wrong, position=f_button.position, text_colour=misc.constants.C_BLACK)
        
        e_button = stimuli.Rectangle(size=(198,50), position=(-100,-300), colour=misc.constants.C_GREY)
        e_text = stimuli.TextLine(text=instructions_feedback_own_answer, position=e_button.position, text_colour=misc.constants.C_BLACK)
        
        r_button = stimuli.Rectangle(size=(198,50), position=(100,-300), colour=misc.constants.C_GREY)
        r_text = stimuli.TextLine(text=instructions_feedback_correct_answer, position=r_button.position, text_colour=misc.constants.C_BLACK)
        
        w_button = stimuli.Rectangle(size=(398,50), position=(0,-352), colour=misc.constants.C_GREY)
        w_text = stimuli.TextLine(text=instructions_next, position=w_button.position, text_colour=misc.constants.C_BLACK)
        
        self.show_grid(c_s1, c_s2, c_s3, c_s4, False)
        
        w_button.plot(self.canvas_stim)
        w_text.plot(self.canvas_stim)
        l_button.plot(self.canvas_stim)
        l_text.plot(self.canvas_stim)
        f_button.plot(self.canvas_stim)
        f_text.plot(self.canvas_stim)
        e_button.plot(self.canvas_stim)
        e_text.plot(self.canvas_stim)
        r_button.plot(self.canvas_stim)
        r_text.plot(self.canvas_stim)
        self.my_control.create_screen(self.canvas_stim, 0, instruction)
        
        return (e_button, r_button, w_button)
    
    def next_task(self, exp, instruction, time):
        button_corr = stimuli.Rectangle(size=(100,50), position=(0,150), colour=misc.constants.C_GREY)
        corr_text = stimuli.TextLine(text=instructions_feedback_correct, position=button_corr.position, text_colour=misc.constants.C_BLACK)
        w_button = stimuli.Rectangle(size=(100,50), position=(0,-300), colour=misc.constants.C_GREY)
        w_text = stimuli.TextLine(text=instructions_next, position=w_button.position, text_colour=misc.constants.C_BLACK)
        w_button.plot(self.canvas_stim)
        w_text.plot(self.canvas_stim)
        button_corr.plot(self.canvas_stim)
        corr_text.plot(self.canvas_stim)
        self.my_control.create_screen(self.canvas_stim, 0, instruction)
        _id, pos, _rt = exp.mouse.wait_press()
        #if w_button.overlapping_with_position(pos):     
        return (self.size1+1, self.size2+1, self.size3+1, self.size_corr1+1, self.size_corr2+1, self.size_corr3+1,  time)
        
        
    def manage_adjustable_stimuli(self, exp, s1,s2,s3,s4, c_s1,c_s2,c_s3,c_s4, bool_instruction,  add_text=[]):
        self.size_corr1 = int(c_s1)-1
        self.size_corr2 = int(c_s2)-1
        self.size_corr3 = int(c_s3)-1
        self.size_corr4 = int(c_s4)-1
        
        time = 0
        if bool_instruction == 2:
            return (self.size1+1, self.size2+1, self.size3+1, self.size4+1, self.size_corr1, self.size_corr2, self.size_corr3, time)
        
        if add_text != []:
            instruction = stimuli.TextBox(text=add_text[0], size = (control.defaults.window_size[0]-100, control.defaults.window_size[1]-100), text_justification = 0, text_colour = misc.constants.C_BLACK)
            self.my_control.create_screen(self.canvas_stim, 0, instruction)
        else:
            instruction = 0
            self.my_control.create_screen(self.canvas_stim, 0, 0)
        
        while True:
            _id, pos, _rt = exp.mouse.wait_press()
            if self.button_ok.overlapping_with_position(pos):
                if (bool_instruction == 1):
                    return (self.size1+1, self.size2+1, self.size3+1, self.size4+1, self.size_corr1, self.size_corr2, self.size_corr3, self.size_corr4, time)
                if (self.size1 == self.size_corr1 and self.size2 == self.size_corr2 and self.size3 == self.size_corr3 and self.size4 == self.size_corr4):
                    return self.next_task(exp, instruction, time)
                
                if not (self.size1 == self.size_corr1 and self.size2 == self.size_corr2 and self.size3 == self.size_corr3 and self.size4 == self.size_corr4):
                    e, c, w = self.give_feedback(self.size_corr1, self.size_corr2, self.size_corr3, instructions_feedback_correct_answer, darkgreen, instruction)
                    while True:
                        _id, new_pos, _rt = exp.mouse.wait_press()
                        if e.overlapping_with_position(new_pos):
                            e, c, w = self.give_feedback(self.size1, self.size2, self.size3, self.size4, instructions_feedback_own_answer, darkred, instruction)
                        if c.overlapping_with_position(new_pos):
                            e, c, w = self.give_feedback(self.size_corr1, self.size_corr2, self.size_corr3, self.size_corr4,  instructions_feedback_correct_answer, darkgreen, instruction)
                        if w.overlapping_with_position(new_pos):
                            return (self.size1+1, self.size2+1, self.size3+1, self.size_corr1+1, self.size_corr2+1, self.size_corr3+1, self.size_corr4+1,  time)                       
            
            
            bool_change = False
            
            i=0
            for cell in self.cells1:
                if cell.overlapping_with_position(pos):
                    time = time+_rt
                    self.size1 = i
                    bool_change = True
                i=i+1      
            
            i=0
            for cell in self.cells2:            
                if cell.overlapping_with_position(pos):
                    time = time+_rt
                    #size2 = self.update_size(i)
                    self.size2 = i
                    bool_change = True
                i=i+1
                    
                    
            i=0
            for cell in self.cells3:
                if cell.overlapping_with_position(pos):
                    time = time+_rt
                    self.size3 = i
                    bool_change = True
                i=i+1
                    
            i=0
            for cell in self.cells4:
                if cell.overlapping_with_position(pos):
                    time = time+_rt
                    self.size4 = i
                    bool_change = True
                i=i+1
                
            if bool_change:
                self.show_grid(self.size1, self.size2, self.size3, self.size4)
                if add_text != []:
                    instruction = stimuli.TextBox(text=add_text[0], size = (control.defaults.window_size[0]-100, control.defaults.window_size[1]-100), text_justification = 0, text_colour = misc.constants.C_BLACK)
                    self.my_control.create_screen(self.canvas_stim, 0, instruction)
                else:
                    instruction = 0
                    self.my_control.create_screen(self.canvas_stim, 0, 0)
                bool_change = True 
                
                  

    def init_cues_withControl(self, my_control, s1, s2, s3, pal_set, pos_set, init = True, bool_ok = False):
        self.my_control = my_control
        arr = (s1,s2,s3)
        self.size1 = int(arr[pos_set[0]])-1
        self.size2 = int(arr[pos_set[1]])-1
        self.size3 = int(arr[pos_set[2]])-1
        self.size4 = None
        self.size_corr1 = int(s1)-1
        self.size_corr2 = int(s2)-1
        self.size_corr3 = int(s3)-1
        self.size_corr4 = None
        self.create_grid_3(pal_set, pos_set)
        self.show_grid(self.size1, self.size2, self.size3, bool_ok)
            
    def init_cues(self, s1, s2, s3, pal_set, pos_set, init = True, bool_ok = False):
        arr = (s1,s2,s3)
        self.size1 = int(arr[pos_set[0]])-1
        self.size2 = int(arr[pos_set[1]])-1
        self.size3 = int(arr[pos_set[2]])-1
        self.size4 = None
        self.size_corr1 = int(s1)-1
        self.size_corr2 = int(s2)-1
        self.size_corr3 = int(s3)-1
        self.size_corr4 = None
        self.create_grid_3(pal_set, pos_set)
        self.show_grid(self.size1, self.size2, self.size3, bool_ok)
