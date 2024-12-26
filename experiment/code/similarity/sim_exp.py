# -*- coding: utf-8 -*-
# Similarity Experiment
# 1. Import necessary libraries
import expyriment
import csv
import os
#from typing import List
import numpy
import numpy.random
#import turtle
import random
import itertools
from stimulus import AdjustableStimulus
from win32api import GetSystemMetrics # Windows: install pywin32
# new packages
import math
import webbrowser
#import subprocess

# 2. Create Exemplars
# 2.1. Generate new class "Exemplar" consisting of name, values, cirt, and pos
class Exemplar:
    def __init__(self, name, values, crit, type, pos=-1):
        self.name = name
        self.values = values
        self.crit = crit
        self.type = type
        self.pos = pos

# 2.2. Define function that reads csv-file containing the exemplars
def read_csv(file, cue_num):
    exemplars = []
    i = 0
    cue_pos_vec=range(cue_num)
    #random.shuffle(cue_pos_vec)
    with open(file, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        for row in reader:
            r = row[0].split(";")
            name = r[0]
            a = (r[1], r[2], r[3])
            if i>0:
                values = (a[cue_pos_vec[0]], a[cue_pos_vec[1]], a[cue_pos_vec[2]])
                crit = r[4]
                type = r[5]
                ex = Exemplar(name, values, -1, type)
                exemplars.append(ex)
            i = i+1
    return exemplars

# 2.3. Read exemplars for experiment
#import rpy2.robjects as robjects
#robjects.r.source("sample.R")
#subprocess.check_output(['R-Datei', "make_stimuli.R"])

exemplars_familiarization_left = read_csv(os.path.join(os.path.dirname(__file__), 'familiarization.stimuli.left.csv'), 3)  # type: List[Exemplar]
exemplars_familiarization_right = read_csv(os.path.join(os.path.dirname(__file__), 'familiarization.stimuli.right.csv'), 3)  # type: List[Exemplar]
random_order_familiarization = numpy.array(range(len(exemplars_familiarization_left)))
random_order_familiarization = random.sample(random_order_familiarization, len(random_order_familiarization))

exemplars_test_1 = read_csv(os.path.join(os.path.dirname(__file__), 'stimuli.1.csv'), 3)  # type: List[Exemplar]
exemplars_test_2 = read_csv(os.path.join(os.path.dirname(__file__), 'stimuli.2.csv'), 3)  # type: List[Exemplar]

no_i_and_v = list((random.choice((1, 2)), random.choice((3, 4))))
print no_i_and_v

# 2.4. Randomization of colour and shape
sets = list(itertools.permutations(range(3)))
pal_nr = random.sample(range(len(sets)), 1)[0]
pos_nr = random.sample(range(len(sets)), 1)[0]
pal_set = sets[pal_nr]
pos_set = sets[pos_nr]
pal_set_var = ''.join(str(e) for e in pal_set)
pos_set_var = ''.join(str(e) for e in pos_set)

print(pal_set_var, pos_set_var)

# 3. Experiment Structure
# Screen Resolution
x = GetSystemMetrics(0)
y = GetSystemMetrics(1)

# 3.1. Experiment initialization
exp = expyriment.design.Experiment(name="sim_exp")
expyriment.io.defaults.outputfile_time_stamp = True # omits time of experiment conduction from file name in data folder
expyriment.control.defaults.window_size = (1500, y-100)
expyriment.control.set_develop_mode(False) # set develop mode so one does not have to render subject id, cool!
expyriment.control.initialize(exp)

exp.mouse.show_cursor()

# 2.5. Preload all possible stimuli
# 3.2. Preload all possible stimuli
stimuli_right = {}
stimuli_left = {}
for i in numpy.array(range(1, 6)):
    for j in numpy.array(range(1, 6)):
        for k in numpy.array(range(1, 6)):
            stim_right = AdjustableStimulus(side="right")
            stim_right.init_cues(i, j, k, pal_set, pos_set, True, False)
            stimuli_right[str(i) + str(j) + str(k)] = stim_right

            stim_left = AdjustableStimulus(side="left")
            stim_left.init_cues(i, j, k, pal_set, pos_set, True, False)
            stimuli_left[str(i) + str(j) + str(k)] = stim_left

# 3.2. Instructions
text_font="Century Gothic"

# 3.2.1. General instructions
general_instructions_page_1 = u"Herzlich Willkommen zu dieser Studie, in der wir untersuchen möchten, wie Menschen Ähnlichkeiten wahrnehmen. \n\n" \
                              u"Bitte lesen Sie die Instruktionen aufmerksam durch. Sollten die Instruktionen unklar sein oder sollte das Experiment nicht richtig funktionieren, geben Sie bitte umgehend dem Studienleiter Bescheid."
general_instructions_page_1 = expyriment.stimuli.TextScreen("", general_instructions_page_1, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

general_instructions_page_2 = u"In dieser Studie geht es darum die Ähnlichkeit zwischen zwei Produkten zu beurteilen. \n\n" \
                              u"Jedes Produkt setzt sich aus drei Zutaten zusammen, welche in unterschiedlicher Menge vorhanden sein können. " \
                              u"Die Produkte unterscheiden sich darin, in welcher Menge die drei Zutaten vorhanden sind. \n\n " \
                              u"Im Folgenden sehen sie zwei Beispielprodukte. Bitte schauen Sie sich die Beispiele in Ruhe an, da es wichtig für die Studie ist, dass sich sich mit den Produkten vertraut machen."
general_instructions_page_2 = expyriment.stimuli.TextScreen("", general_instructions_page_2, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

general_instructions_page_3 = u"Hier sehen Sie das erste Beispielprodukt. Jede Zutat ist schematisch durch einen Balken dargestellt. " \
                              u"Die Anzahl farbiger Einheiten innerhalb der Balken gibt an, in welcher Menge die jeweilige Zutat vorhanden ist. \n\n" \
                              u"Das hier dargestellte Produkt hat von jeder Zutat die Minimalmenge 1. \n" \
                              u"Es ist unmöglich ein Produkt mit weniger als einer Mengeneinheit einer Zutat herzustellen. "
general_instructions_page_3 = expyriment.stimuli.TextScreen("", general_instructions_page_3, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

general_instructions_page_4 = u"Hier sehen Sie das zweite Beispielprodukt. Jede Zutat ist schematisch durch einen Balken dargestellt. " \
                              u"Die Anzahl farbiger Einheiten innerhalb der Balken gibt an, in welcher Menge die jeweilige Zutat vorhanden ist. \n\n" \
                              u"Das hier dargestellte Produkt hat von jeder Zutat die Maximalmenge 5. \n" \
                              u"Es ist unmöglich ein Produkt mit mehr als fünf Mengeneinheiten einer Zutat herzustellen. "
general_instructions_page_4 = expyriment.stimuli.TextScreen("", general_instructions_page_4, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

general_instructions_page_5 = u"In jedem Durchgang werden Ihnen zwei zufällig ausgewählte Produkte gezeigt. \n\n" \
                              u"Ihre Aufgabe ist es, die Ähnlichkeit der zwei Produkte zu beurteilen. Ihr Urteil können Sie eingeben, indem Sie auf den Balken klicken."
general_instructions_page_5 = expyriment.stimuli.TextScreen("", general_instructions_page_5, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

general_instructions_similar = u"Je weiter rechts Sie auf den Balken klicken, desto ähnlicher empfinden Sie die zwei Produkte."
general_instructions_similar = expyriment.stimuli.TextScreen("", general_instructions_similar, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

general_instructions_dissimilar = u"Je weiter links Sie auf den Balken klicken, desto unähnlicher empfinden Sie die zwei Produkte."
general_instructions_dissimilar = expyriment.stimuli.TextScreen("", general_instructions_dissimilar, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

general_instructions_identical = u"Wenn Sie zwei identische Produkte sehen, klicken Sie bitte ganz rechts auf den Balken. " \
                                 u"Bei allen Paaren von Produkten, die nicht identisch sind, zählt Ihre subjektive Einschätzung und es gibt somit kein richtig oder falsch."
general_instructions_identical = expyriment.stimuli.TextScreen("", general_instructions_identical, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

general_instructions_page_6 = u"Hier sehen Sie zwei identische Produkte. Klicken Sie bitte ganz rechts in den roten Kasten auf dem Balken."
general_instructions_page_6 = expyriment.stimuli.TextScreen("", general_instructions_page_6, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

# 3.2.1.2. Time Pressure-Bedingung (time-pressure_first vs. time-pressure_last)
general_instructions_page_7 = u"Im Laufe der Studie werden Sie die Ähnlichkeit zwischen zwei Produkten in drei Phasen beurteilen. " \
                              u"Die Phasen des Experiments unterscheiden sich in den Produktpaaren und darin, wieviel Zeit Sie zur Beurteilung der Ähnlichkeit haben."
general_instructions_page_7 = expyriment.stimuli.TextScreen("", general_instructions_page_7, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

# 3.2.1.1. Familiarization instructions
familiarization_instructions = u"Sie können nun mit der ersten Phase beginnen. \n\n" \
                               u"In dieser Phase können Sie die Produkte in jedem Durchgang so lange anschauen, wie Sie möchten. " \
                               u"Jeder Durchgang beginnt damit, dass zwei Produkte eingeblendet werden. Sie können sich dann die Produkte in Ruhe anschauen und überlegen wie ähnlich sie sich sehen. \n" \
                               u"Wenn Sie zu einem Urteil über die Ähnlichkeit gelangt sind, klicken Sie bitte auf den Weiter-Knopf. Dann erscheint, der  Balken, auf dem Sie Ihre Antwort eingeben können. " \
                               u"Bitte beachten Sie, dass Sie Ihre Antwort innerhalb von 3 Sekunden auf dem Balken eingeben müssen, damit der Durchgang zählt. \n\n" \
                               u"Denken Sie daran, je weiter rechts Sie auf den Balken klicken, desto ähnlicher empfinden Sie die zwei Produkte. " \
                               u"Je weiter links Sie auf den Balken klicken, desto unähnlicher empfinden Sie die zwei Produkte. \n\n" \
                               u"Nachdem Sie Ihre Antwort eingegeben haben, können Sie durch Klicken auf den Weiter-Knopf den nächsten Durchgang starten. "
familiarization_instructions = expyriment.stimuli.TextScreen("", familiarization_instructions, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

false_click_familiarization = u"Klicken Sie bitte ganz rechts in den roten Kasten auf dem Balken."
false_click_familiarization = expyriment.stimuli.TextScreen("", false_click_familiarization, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), text_bold=True, size = (0.75*x, y/2), text_font=text_font)

false_click = u"Klicken Sie bitte auf den Balken um Ihre Antwort zu geben."
false_click= expyriment.stimuli.TextScreen("", false_click, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), text_bold=True, size = (0.75*x, y/2), text_font=text_font)

# 3.2. Testphasen
tp_fam = u"Sie haben zwischen %d und %d Sekunden Zeit sich die Produkte anzuschauen und zwischen %d und %d Sekunden Zeit Ihre Antwort einzugeben. \n\n" \
         u"Klicken Sie auf Weiter um mit den Übungsdurchgängen zu beginnen."

tp_fam_post = u"Sie haben die Übungsdurchgänge nun absolviert und können jetzt mit den eigentlichen Durchgängen beginnen."
tp_fam_post = expyriment.stimuli.TextScreen("", tp_fam_post, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

# 3.2.1. Time pressure first
tp_first_test1 = u"Nun beginnt die zweite Phase des Experiments. \n\n" \
                 u"Diese Phase unterscheidet sich von der vorherigen Phase darin, dass Sie nicht nur begrenzt Zeit haben, sich die Produkte anzuschauen, sondern auch weniger Zeit haben Ihre Antwort einzugeben. Sonst bleibt alles gleich. \n\n" \
                 u"Bitte beachten Sie, dass Sie sowohl rechtzeitig auf den Weiter-Knopf drücken müssen um zum Antwortbalken zu gelangen als auch rechtzeitig Ihre Antwort auf dem Balken angeben müssen, damit der Durchgang zählt. \n\n" \
                 u"Es ist wichtig, dass Sie das Zeitlimit nicht überschreiten! \n\n" \
                u"Im Folgenden können Sie in einigen Übungsdurchgängen lernen, wie viel Zeit Sie haben, und sich an das Zeitlimit gewöhnen."
tp_first_test1 = expyriment.stimuli.TextScreen("", tp_first_test1, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)


tp_first_test2 = u"Nun beginnt die letzte Phase des Experiments. \n\n" \
                u"In dieser Phase haben Sie soviel Zeit, wie sie wollen, um sich die Produkte anzuschauen und deren Ähnlichkeit zu beurteilen. \n\n" \
                u"Wie gehabt, gelangen Sie durch Klicken auf den Weiter-Knopf zum Balken, auf dem Sie Ihre Antwort eingeben können. \n" \
                u"Sie haben wieder 3 Sekunden um Ihre Antwort auf dem Balken einzugeben. "
tp_first_test2 = expyriment.stimuli.TextScreen("", tp_first_test2, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

# 3.2.2. Time pressure last
tp_last_test1 = u"Nun beginnt die zweite Phase des Experiments. \n\n" \
                u"In dieser Phase haben Sie soviel Zeit, wie sie wollen, um sich die Produkte anzuschauen und deren Ähnlichkeit zu beurteilen. \n\n" \
                u"Wie gehabt, gelangen Sie durch Klicken auf den Weiter-Knopf zum Balken, auf dem Sie Ihre Antwort eingeben können. \n" \
                u"Sie haben wieder 3 Sekunden um Ihre Antwort auf dem Balken einzugeben. "
tp_last_test1 = expyriment.stimuli.TextScreen("", tp_last_test1, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

tp_last_test2 = u"Nun beginnt die letzte Phase des Experiments. \n\n" \
                 u"Diese Phase unterscheidet sich von der vorherigen Phase darin, dass Sie nicht nur begrenzt Zeit haben, sich die Produkte anzuschauen, sondern auch weniger Zeit haben Ihre Antwort einzugeben. Sonst bleibt alles gleich. \n\n" \
                 u"Bitte beachten Sie, dass Sie sowohl rechtzeitig auf den Weiter-Knopf drücken müssen um zum Antwortbalken zu gelangen als auch rechtzeitig Ihre Antwort auf dem Balken angeben müssen, damit der Durchgang zählt. \n\n" \
                 u"Es ist wichtig, dass Sie das Zeitlimit nicht überschreiten! \n\n" \
                u"Im Folgenden können Sie in einigen Übungsdurchgängen lernen, wie viel Zeit Sie haben, und sich an das Zeitlimit gewöhnen."
tp_last_test2 = expyriment.stimuli.TextScreen("", tp_last_test2, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

# 3.2.3. Post-test instruction
post_test_instructions = u"Sie haben die letzte Phase nun abgeschlossen. \n\n" \
                         u"Durch Drücken des Weiter-Knopfs gelangen Sie zum demographischen Fragebogen. \n\nBitte füllen Sie diesen zum Schluss noch aus."
post_test_instructions = expyriment.stimuli.TextScreen("", post_test_instructions, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y/2), text_font=text_font)

# 3.2.4. Instructions if too slow
too_slow_instructions = u"Sie waren zu langsam! Versuchen Sie das nächste Mal schneller zu antworten!"
too_slow_instructions = expyriment.stimuli.TextScreen("", too_slow_instructions, text_bold=True, text_colour=(0, 0, 0), background_colour=(255, 255, 255), text_size=24, position=(50,50), size = (0.75*x, y/2), text_font=text_font)

# 3.3. Conditions
bool_timepressure_first = False
contemplation_time_vec = []
reaction_time_vec = []

# 5. Define writer
writer = None

# 6. Response slider
# Global line to check for position overlap
line = expyriment.stimuli.Line(start_point=(-500,0), end_point=(500,0), line_width=60)
line.preload()

# Slider with color gradient
steps = 20
lines = {}
for i in range(0, steps):
    lines[str(i)] = expyriment.stimuli.Line(start_point=(i*1000/steps - 500, 0), end_point=(i*1000/steps - (500 - 1000/steps), 0), line_width=30, colour = (-i*200/steps + 200, -i*200/steps + 200, -i*200/steps + 200))

# Slider ticks
left_line = expyriment.stimuli.Line(start_point = (-501, -30), end_point=(-501,30), line_width = 2, colour = (0, 0, 0))
middle_line = expyriment.stimuli.Line(start_point = (0, -30), end_point=(0,30), line_width = 2, colour = (0, 0, 0))
right_line = expyriment.stimuli.Line(start_point = (501, -30), end_point=(501,30), line_width = 2, colour = (0, 0, 0))

# Slider anchor labels
identical = expyriment.stimuli.TextBox("identisch", size = (0.1*x, 50), position = (line.end_point[0] + 0.06*x, -10), text_justification=0,
                                       text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), text_font = text_font)
identical.preload()

different = expyriment.stimuli.TextBox("komplett\nverschieden", size = (0.1*x, 80), position = (line.start_point[0] - 0.06*x, -10), text_justification=2,
                                       text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), text_font = text_font)
different.preload()

# Global line for instructions
line_instructions = expyriment.stimuli.Line(start_point=(-500,-250), end_point=(500,-250), line_width=30)
line_instructions.preload()

lines_instructions = {}
for i in range(0, steps):
    lines_instructions[str(i)] = expyriment.stimuli.Line(start_point=(i*1000/steps - 500, -250), end_point=(i*1000/steps - (500 - 1000/steps), -250), line_width=30, colour = (-i*200/steps + 200, -i*200/steps + 200, -i*200/steps + 200))

left_line_instructions = expyriment.stimuli.Line(start_point = (-501, -280), end_point=(-501,-220), line_width = 2, colour = (0, 0, 0))
middle_line_instructions = expyriment.stimuli.Line(start_point = (0, -280), end_point=(0,-220), line_width = 3, colour = (0, 0, 0))
right_line_instructions = expyriment.stimuli.Line(start_point = (501, -280), end_point=(501,-220), line_width = 2, colour = (0, 0, 0))

identical_instructions = expyriment.stimuli.TextBox("identisch", size = (0.1*x, 50), position = (line.end_point[0] + 0.06*x, -260), text_justification=0,
                                       text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), text_font = text_font)
identical_instructions.preload()

different_instructions = expyriment.stimuli.TextBox("komplett\nverschieden", size = (0.1*x, 80), position = (line.start_point[0] - 0.06*x, -260), text_justification=2,
                                       text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), text_font = text_font)
different_instructions.preload()

# 7. Create arrows and rectangles for instruction
arrow_right = expyriment.stimuli.Shape(position=(350, -100), colour=(255, 0 ,0), line_width=30)
arrow_right.add_vertices([(-50, 50), (50, -50), (-150, 0), (150, 0), (-50, -50)])

arrow_left = expyriment.stimuli.Shape(position=(-350, -100), colour=(255, 0 ,0), line_width=30)
arrow_left.add_vertices([(50, 50), (-50, -50), (150, 0), (-150, 0), (50, -50)])

rectangle = expyriment.stimuli.Rectangle(position=(475, -250), colour=(255, 0 ,0), size = (50, 30), line_width=5)
rectangle.preload()

rectangle_2 = expyriment.stimuli.Rectangle(position=(475, -250), colour=(255, 0 ,0), size = (50, 60))
rectangle_2.preload()

# 8. Weiter Button
def wait_for_next(canvas, crit=-1, rt=-1, next_bool=True, pos = (0, (-y/2+100)), contemplation_time = -1, instructions = False):
    button_next = expyriment.stimuli.Rectangle(size=(200, 50), colour= expyriment.misc.constants.C_GREY, line_width=0, position=pos)
    text_next = expyriment.stimuli.TextLine("Weiter", pos, text_colour = expyriment.misc.constants.C_BLACK, text_font = text_font)
    button_next.plot(canvas)
    text_next.plot(canvas)
    canvas.present()
    if instructions == True:
        exp.clock.wait(1000)	
    if contemplation_time == -1:
        id, pos_press, rt = exp.mouse.wait_press()
    else:
        id, pos_press, rt = exp.mouse.wait_press(duration = contemplation_time)
    if button_next.overlapping_with_position(pos_press):
        return rt
    else:
        if contemplation_time == -1:
            return wait_for_next(canvas, crit, rt, next_bool, pos, contemplation_time)
        else:
            return rt

# 9. Fixation cross
fixation_cross = expyriment.stimuli.Shape(position = (0, 0), colour = (0, 0, 0), line_width=2)
fixation_cross.add_vertices([(25, 0), (-50, 0), (25, 0), (0, 25), (0, -50), (0, 25)])


########################################################################################################################
##################################################  Start Experiment ###################################################
########################################################################################################################
expyriment.control.start()

# Write data into participant specific result file
with open(str("data/results_" + str(exp.subject) + ".csv"), "wb") as csvfile:
    writer = csv.writer(csvfile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    # Header
    row = ["subj_id", "trial", "block", "stim_left", "stim_right", "type", "response",
           "contemplation_time", "reaction_time", "contemplation_time_pressure", "reaction_time_pressure",
           "color_presentation_order", "shape_presentation_order", "time_pressure_first"]
    writer.writerow(row)

    # Condition assignment: timepressure_first (True or False)
    if (exp.subject % 2) == 1:
        bool_timepressure_first = True

    # Define stimuli as adjustable stimuli and canvas as blank screen
    canvas = expyriment.stimuli.Canvas(size=(x, y),colour=expyriment.misc.constants.C_WHITE)
    stim_middle = AdjustableStimulus(side="middle")
    stim_right = AdjustableStimulus(side="right")
    stim_left = AdjustableStimulus(side="left")


    ####################################################################################################################
    ###############################################  Global Instructions ###############################################
    ####################################################################################################################
    general_instructions_page_1.plot(canvas)
    wait_for_next(canvas, instructions = True)

    general_instructions_page_2.plot(canvas)
    wait_for_next(canvas, instructions = True)

    general_instructions_page_3.plot(canvas)
    stim_middle.init_cues(1, 1, 1, pal_set, pos_set, True, False)
    stim_middle.canvas_stim.plot(canvas)
    wait_for_next(canvas, instructions = True)

    general_instructions_page_4.plot(canvas)
    stim_middle.init_cues(5, 5, 5, pal_set, pos_set, True, False)
    stim_middle.canvas_stim.plot(canvas)
    wait_for_next(canvas, instructions = True)

    general_instructions_page_5.plot(canvas)
    left_line_instructions.plot(canvas)
    middle_line_instructions.plot(canvas)
    right_line_instructions.plot(canvas)
    [lines_instructions[str(j)].plot(canvas) for j in range(0, len(lines_instructions))]
    identical_instructions.plot(canvas)
    different_instructions.plot(canvas)
    wait_for_next(canvas, instructions = True)

    general_instructions_similar.plot(canvas)
    left_line_instructions.plot(canvas)
    middle_line_instructions.plot(canvas)
    right_line_instructions.plot(canvas)
    [lines_instructions[str(j)].plot(canvas) for j in range(0, len(lines_instructions))]
    identical_instructions.plot(canvas)
    different_instructions.plot(canvas)
    arrow_right.plot(canvas)
    wait_for_next(canvas, instructions = True)

    general_instructions_dissimilar.plot(canvas)
    left_line_instructions.plot(canvas)
    middle_line_instructions.plot(canvas)
    right_line_instructions.plot(canvas)
    [lines_instructions[str(j)].plot(canvas) for j in range(0, len(lines_instructions))]
    identical_instructions.plot(canvas)
    different_instructions.plot(canvas)
    arrow_left.plot(canvas)
    wait_for_next(canvas, instructions = True)

    canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
    general_instructions_identical.plot(canvas)
    wait_for_next(canvas, instructions = True)

    valid_pos = False
    while valid_pos == False:
        # Show two identical stimuli and a line with a rectangle showing where people need to click in that case
        canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
        general_instructions_page_6.plot(canvas)
        stim = stimuli_left["423"]
        stim.canvas_stim.plot(canvas)
        stim = stimuli_right["423"]
        stim.canvas_stim.plot(canvas)
        left_line_instructions.plot(canvas)
        middle_line_instructions.plot(canvas)
        right_line_instructions.plot(canvas)
        [lines_instructions[str(j)].plot(canvas) for j in range(0, len(lines_instructions))]
        identical_instructions.plot(canvas)
        different_instructions.plot(canvas)
        rectangle.plot(canvas)
        canvas.present()

        # Response
        click = exp.mouse.wait_press(0)
        pos = click[1]

        # Check whether click overlaps with predefined rectangle
        if rectangle_2.overlapping_with_position(pos) == False:
            canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
            false_click_familiarization.plot(canvas)
            canvas.present()
            exp.clock.wait(1500)
        else:
            valid_pos = True
    # If click was in rectangle present wait_for_next
    wait_for_next(canvas)
    canvas.present(clear = False)

    general_instructions_page_7.plot(canvas)
    wait_for_next(canvas, instructions = True)


    ####################################################################################################################
    ##############################################  Familiarization Phase ##############################################
    ####################################################################################################################

    # Instructions
    familiarization_instructions.plot(canvas)
    wait_for_next(canvas, instructions = True)

    for i in random_order_familiarization:
        canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
        fixation_cross.plot(canvas)
        canvas.present()
        exp.clock.wait(500)

        # Contemplation subphase
        canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
        id_left = exemplars_familiarization_left[i].name
        stim = stimuli_left[id_left]
        stim.canvas_stim.plot(canvas)

        id_right = exemplars_familiarization_right[i].name
        stim = stimuli_right[id_right]
        stim.canvas_stim.plot(canvas)

        rt = wait_for_next(canvas)

        contemplation_time = rt
        contemplation_time_vec.append(contemplation_time)

        # Response subphase
        valid_pos = False
        rt = 0
        while valid_pos == False:
            # Display slider
            canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
            left_line.plot(canvas)
            middle_line.plot(canvas)
            right_line.plot(canvas)
            [lines[str(j)].plot(canvas) for j in range(0, len(lines))]
            identical.plot(canvas)
            different.plot(canvas)
            canvas.present()

            # Response
            click = exp.mouse.wait_press(0, duration=3000-rt)

            # Check whether no click has been made
            if click[2] == None:
                canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                too_slow_instructions.plot(canvas)
                canvas.present()
                pos = ["NA"]
                rt = 3000
                valid_pos = True
            else:
                pos = click[1]
                rt = rt + click[2]

                # Check whether click was made on the slider
                if line.overlapping_with_position(pos):
                    marker = expyriment.stimuli.Line(
                        start_point=(pos[0], line.start_point[1] + line.line_width/2 * -0.75),
                        end_point=(pos[0], line.start_point[1] + line.line_width/2 * 0.75),
                        line_width=12,
                        colour=(255, 0, 0))
                    marker.plot(canvas)
                    canvas.present(clear=False)
                    valid_pos = True
                else:
                    canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                    false_click.plot(canvas)
                    canvas.present()
                    exp.clock.wait(1500)

        wait_for_next(canvas)
        canvas.present(clear=False)

        # Eintrag für einen trial
        row = [exp.subject, i, "familiarization", id_left, id_right, "fam", pos[0],
               contemplation_time, rt, "NA", "NA",
               pal_set_var, pos_set_var, bool_timepressure_first]
        writer.writerow(row)

        reaction_time_vec.append(rt)

    # Determine Time Pressure
    time = numpy.median(contemplation_time_vec[-15:]) * .15
    reaction_time = numpy.percentile(reaction_time_vec[-15:], 90)

    ####################################################################################################################
    ################################################### Test Phase 1 ###################################################
    ####################################################################################################################

    # Familiarization with time pressure
    canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
    if bool_timepressure_first:
        tp_first_test1.plot(canvas)
        wait_for_next(canvas, instructions = True)

        tp_fam = tp_fam % (math.floor(time / 1000.0), math.ceil(time / 1000.0), math.floor(reaction_time / 1000.0),
                           math.ceil(reaction_time / 1000.0))
        tp_fam = expyriment.stimuli.TextScreen("", tp_fam, text_size=24, background_colour=(255, 255, 255),
                                               text_colour=(0, 0, 0), size=(0.75*x, y), text_font=text_font)
        tp_fam.plot(canvas)
        wait_for_next(canvas, instructions = True)

        for i in random_order_familiarization:
            canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
            fixation_cross.plot(canvas)
            canvas.present()
            exp.clock.wait(500)

            # Contemplation subphase
            canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
            id_left = exemplars_familiarization_left[i].name
            stim = stimuli_left[id_left]
            stim.canvas_stim.plot(canvas)

            id_right = exemplars_familiarization_right[i].name
            stim = stimuli_right[id_right]
            stim.canvas_stim.plot(canvas)

            contemplation_time = wait_for_next(canvas, contemplation_time = time)

            # Check if people did not press wait_for_next
            if contemplation_time is None:
                canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                too_slow_instructions.plot(canvas)
                canvas.present()
                rt = "NA"
                pos = ["NA"]
                valid_pos = True

                # canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                # rt = wait_for_next(canvas, contemplation_time = reaction_time)
                #
                # if rt is None:
                #     canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                #     too_slow_instructions.plot(canvas)
                #     canvas.present()
                #     rt = "NA"
                #     pos = ["NA"]
                #     valid_pos = True
                # else:
                #     canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                #     left_line.plot(canvas)
                #     middle_line.plot(canvas)
                #     right_line.plot(canvas)
                #     [lines[str(j)].plot(canvas) for j in range(0, len(lines))]
                #     identical.plot(canvas)
                #     different.plot(canvas)
                #     canvas.present()
                #     valid_pos = False
                #     rt = 0 + rt
            else:
                canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                left_line.plot(canvas)
                middle_line.plot(canvas)
                right_line.plot(canvas)
                [lines[str(j)].plot(canvas) for j in range(0, len(lines))]
                identical.plot(canvas)
                different.plot(canvas)
                canvas.present()
                valid_pos = False
                rt = 0

            while valid_pos == False and rt < reaction_time:
                click = exp.mouse.wait_press(0, duration=reaction_time - rt)
                pos = click[1]
                valid_pos = line.overlapping_with_position(pos)
                if click[2] != None:
                    rt = rt + click[2]
                else:
                    canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                    too_slow_instructions.plot(canvas)
                    canvas.present()
                    rt = "NA"
                    pos = ["NA"]
            if pos != ["NA"]:
                marker = expyriment.stimuli.Line(
                    start_point=(pos[0], line.start_point[1] + line.line_width/2 * -0.75),
                    end_point=(pos[0], line.start_point[1] + line.line_width/2 * 0.75),
                    line_width=12,
                    colour=(255, 0, 0))
                marker.plot(canvas)
                canvas.present(clear=False)

            wait_for_next(canvas)
            canvas.present(clear = False)

            # Eintrag für einen trial
            row = [exp.subject, i, "familiarization_time_pressure", id_left, id_right, "fam", pos[0],
                   contemplation_time, rt, time, reaction_time,
                   pal_set_var, pos_set_var, bool_timepressure_first]
            writer.writerow(row)

        canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
        tp_fam_post.plot(canvas)
        wait_for_next(canvas, instructions = True)
    else:
        tp_last_test1.plot(canvas)
        wait_for_next(canvas, instructions = True)

    for g in [1, 2, 3, 4]:
        if g in no_i_and_v: # omits type I and V pairs
            exemplars_test_1_temp = exemplars_test_1[:-9]
            exemplars_test_2_temp = exemplars_test_2[:-9]
        else:
            exemplars_test_1_temp = exemplars_test_1
            exemplars_test_2_temp = exemplars_test_2
        print len(exemplars_test_1_temp)

        random_order_test = numpy.array(range(len(exemplars_test_1_temp)))
        random_order_test = random.sample(random_order_test, len(random_order_test))
        for i in random_order_test:
            canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
            fixation_cross.plot(canvas)
            canvas.present()
            exp.clock.wait(500)

            # Contemplation subphase
            canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)

            if (exp.subject % 2) == 1:
                id_left = exemplars_test_1_temp[i].name
                id_right = exemplars_test_2_temp[i].name
            else:
                id_left = exemplars_test_2_temp[i].name
                id_right = exemplars_test_1_temp[i].name
            stim = stimuli_left[id_left]
            stim.canvas_stim.plot(canvas)

            stim = stimuli_right[id_right]
            stim.canvas_stim.plot(canvas)

            if bool_timepressure_first:
                contemplation_time = wait_for_next(canvas, contemplation_time = time)

                if contemplation_time is None:
                    canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                    too_slow_instructions.plot(canvas)
                    canvas.present()
                    rt = "NA"
                    pos = ["NA"]
                    valid_pos = True

                    # canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                    # rt = wait_for_next(canvas, contemplation_time = reaction_time)
                    # if rt is None:
                    #     canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                    #     too_slow_instructions.plot(canvas)
                    #     canvas.present()
                    #     rt = "NA"
                    #     pos = ["NA"]
                    #     valid_pos = True
                    # else:
                    #     canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                    #     left_line.plot(canvas)
                    #     middle_line.plot(canvas)
                    #     right_line.plot(canvas)
                    #     [lines[str(j)].plot(canvas) for j in range(0, len(lines))]
                    #     identical.plot(canvas)
                    #     different.plot(canvas)
                    #     canvas.present()
                    #     valid_pos = False
                    #     rt = 0 + rt
                else:
                    canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                    left_line.plot(canvas)
                    middle_line.plot(canvas)
                    right_line.plot(canvas)
                    [lines[str(j)].plot(canvas) for j in range(0, len(lines))]
                    identical.plot(canvas)
                    different.plot(canvas)
                    canvas.present()
                    valid_pos = False
                    rt = 0
            else:
                contemplation_time = wait_for_next(canvas)

                canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                left_line.plot(canvas)
                middle_line.plot(canvas)
                right_line.plot(canvas)
                [lines[str(j)].plot(canvas) for j in range(0, len(lines))]
                identical.plot(canvas)
                different.plot(canvas)
                canvas.present()
                valid_pos = False
                rt = 0

            if bool_timepressure_first:
                while valid_pos == False and rt < reaction_time:
                    click = exp.mouse.wait_press(0, duration = reaction_time-rt)
                    pos = click[1]
                    valid_pos = line.overlapping_with_position(pos)
                    if click[2] != None:
                        rt = rt + click[2]
                    else:
                        canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                        too_slow_instructions.plot(canvas)
                        canvas.present()
                        rt = "NA"
                        pos = ["NA"]
                if pos != ["NA"]:
                    marker = expyriment.stimuli.Line(
                        start_point=(pos[0], line.start_point[1] + line.line_width/2 * -0.75),
                        end_point=(pos[0], line.start_point[1] + line.line_width/2 * 0.75),
                        line_width=12,
                        colour=(255, 0, 0))
                    marker.plot(canvas)
                    canvas.present(clear=False)

            else:
                while valid_pos == False:
                    click = exp.mouse.wait_press(0)
                    pos = click[1]
                    rt = rt + click[2]
                    valid_pos = line.overlapping_with_position(pos)
                marker = expyriment.stimuli.Line(
                    start_point=(pos[0], line.start_point[1] + line.line_width/2 * -0.75),
                    end_point=(pos[0], line.start_point[1] + line.line_width/2 * 0.75),
                    line_width=12,
                    colour=(255, 0, 0))
                marker.plot(canvas)
                canvas.present(clear=False)

            wait_for_next(canvas)
            canvas.present(clear = False)

            # Eintrag für einen trial
            row = [exp.subject, i, "test_1", id_left, id_right, exemplars_test_1_temp[i].type, pos[0],
                   contemplation_time, rt, time, reaction_time,
                   pal_set_var, pos_set_var, bool_timepressure_first]
            writer.writerow(row)

    ####################################################################################################################
    ################################################### Test Phase 2 ###################################################
    ####################################################################################################################

    canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
    if bool_timepressure_first == False:
        tp_last_test2.plot(canvas)
        wait_for_next(canvas, instructions = True)

        tp_fam = tp_fam % (math.floor(time / 1000.0), math.ceil(time / 1000.0), math.floor(reaction_time / 1000.0),
                           math.ceil(reaction_time / 1000.0))
        tp_fam = expyriment.stimuli.TextScreen("", tp_fam, text_size=24, background_colour=(255, 255, 255),
                                               text_colour=(0, 0, 0), size=(0.75 * x, y), text_font=text_font)
        tp_fam.plot(canvas)
        wait_for_next(canvas, instructions = True)

        for i in random_order_familiarization:
            canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
            fixation_cross.plot(canvas)
            canvas.present()
            exp.clock.wait(500)

            # Contemplation subphase
            canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
            id_left = exemplars_familiarization_left[i].name
            stim = stimuli_left[id_left]
            stim.canvas_stim.plot(canvas)

            id_right = exemplars_familiarization_right[i].name
            stim = stimuli_right[id_right]
            stim.canvas_stim.plot(canvas)

            contemplation_time = wait_for_next(canvas, contemplation_time = time)

            if contemplation_time is None:
                canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                too_slow_instructions.plot(canvas)
                canvas.present()
                rt = "NA"
                pos = ["NA"]
                valid_pos = True

                # canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                # rt = wait_for_next(canvas, contemplation_time = reaction_time)
                # if rt is None:
                #     canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                #     too_slow_instructions.plot(canvas)
                #     canvas.present()
                #     rt = "NA"
                #     pos = ["NA"]
                #     valid_pos = True
                # else:
                #     canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                #     left_line.plot(canvas)
                #     middle_line.plot(canvas)
                #     right_line.plot(canvas)
                #     [lines[str(j)].plot(canvas) for j in range(0, len(lines))]
                #     identical.plot(canvas)
                #     different.plot(canvas)
                #     canvas.present()
                #     valid_pos = False
                #     rt = 0 + rt
            else:
                canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                left_line.plot(canvas)
                middle_line.plot(canvas)
                right_line.plot(canvas)
                [lines[str(j)].plot(canvas) for j in range(0, len(lines))]
                identical.plot(canvas)
                different.plot(canvas)
                canvas.present()
                valid_pos = False
                rt = 0

            while valid_pos == False and rt < reaction_time:
                click = exp.mouse.wait_press(0, duration=reaction_time - rt)
                pos = click[1]
                valid_pos = line.overlapping_with_position(pos)
                if click[2] != None:
                    rt = rt + click[2]
                else:
                    canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                    too_slow_instructions.plot(canvas)
                    canvas.present()
                    rt = "NA"
                    pos = ["NA"]
            if pos != ["NA"]:
                marker = expyriment.stimuli.Line(
                    start_point=(pos[0], line.start_point[1] + line.line_width/2 * -0.75),
                    end_point=(pos[0], line.start_point[1] + line.line_width/2 * 0.75),
                    line_width=12,
                    colour=(255, 0, 0))
                marker.plot(canvas)
                canvas.present(clear=False)

            wait_for_next(canvas)
            canvas.present(clear = False)

            # Eintrag für einen trial
            row = [exp.subject, i, "familiarization_time_pressure", id_left, id_right, "fam", pos[0],
                   contemplation_time, rt, time, reaction_time,
                   pal_set_var, pos_set_var, bool_timepressure_first]
            writer.writerow(row)

        canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
        tp_fam_post.plot(canvas)
        wait_for_next(canvas, instructions = True)
    else:
        tp_first_test2.plot(canvas)
        wait_for_next(canvas, instructions = True)

    for g in [1, 2, 3, 4]:
        if g in no_i_and_v: # omits type I and V pairs
            exemplars_test_1_temp = exemplars_test_1[:-9]
            exemplars_test_2_temp = exemplars_test_2[:-9]
        else:
            exemplars_test_1_temp = exemplars_test_1
            exemplars_test_2_temp = exemplars_test_2
        print len(exemplars_test_1_temp)

        random_order_test = numpy.array(range(len(exemplars_test_1_temp)))
        random_order_test = random.sample(random_order_test, len(random_order_test))
        for i in random_order_test:
            canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
            fixation_cross.plot(canvas)
            canvas.present()
            exp.clock.wait(500)

            # Contemplation subphase
            canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)

            if (exp.subject % 2) == 1:
                id_left = exemplars_test_1_temp[i].name
                id_right = exemplars_test_2_temp[i].name
            else:
                id_left = exemplars_test_2_temp[i].name
                id_right = exemplars_test_1_temp[i].name
            stim = stimuli_left[id_left]
            stim.canvas_stim.plot(canvas)

            stim = stimuli_right[id_right]
            stim.canvas_stim.plot(canvas)

            if bool_timepressure_first == False:
                contemplation_time = wait_for_next(canvas, contemplation_time=time)

                if contemplation_time is None:
                    canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                    too_slow_instructions.plot(canvas)
                    canvas.present()
                    rt = "NA"
                    pos = ["NA"]
                    valid_pos = True

                    # canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                    # rt = wait_for_next(canvas, contemplation_time=reaction_time)
                    # if rt is None:
                    #     canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                    #     too_slow_instructions.plot(canvas)
                    #     canvas.present()
                    #     rt = "NA"
                    #     pos = ["NA"]
                    #     valid_pos = True
                    # else:
                    #     canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                    #     left_line.plot(canvas)
                    #     middle_line.plot(canvas)
                    #     right_line.plot(canvas)
                    #     [lines[str(j)].plot(canvas) for j in range(0, len(lines))]
                    #     identical.plot(canvas)
                    #     different.plot(canvas)
                    #     canvas.present()
                    #     valid_pos = False
                    #     rt = 0 + rt
                else:
                    canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                    left_line.plot(canvas)
                    middle_line.plot(canvas)
                    right_line.plot(canvas)
                    [lines[str(j)].plot(canvas) for j in range(0, len(lines))]
                    identical.plot(canvas)
                    different.plot(canvas)
                    canvas.present()
                    valid_pos = False
                    rt = 0
            else:
                contemplation_time = wait_for_next(canvas)

                canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                left_line.plot(canvas)
                middle_line.plot(canvas)
                right_line.plot(canvas)
                [lines[str(j)].plot(canvas) for j in range(0, len(lines))]
                identical.plot(canvas)
                different.plot(canvas)
                canvas.present()
                valid_pos = False
                rt = 0

            if bool_timepressure_first == False:
                while valid_pos == False and rt < reaction_time:
                    click = exp.mouse.wait_press(0, duration = reaction_time-rt)
                    pos = click[1]
                    valid_pos = line.overlapping_with_position(pos)
                    if click[2] != None:
                        rt = rt + click[2]
                    else:
                        canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                        too_slow_instructions.plot(canvas)
                        canvas.present()
                        rt = "NA"
                        pos = ["NA"]
                if pos != ["NA"]:
                    marker = expyriment.stimuli.Line(
                        start_point=(pos[0], line.start_point[1] + line.line_width/2 * -0.75),
                        end_point=(pos[0], line.start_point[1] + line.line_width/2 * 0.75),
                        line_width=12,
                        colour=(255, 0, 0))
                    marker.plot(canvas)
                    canvas.present(clear=False)

            else:
                while valid_pos == False:
                    click = exp.mouse.wait_press(0)
                    pos = click[1]
                    rt = rt + click[2]
                    valid_pos = line.overlapping_with_position(pos)
                marker = expyriment.stimuli.Line(
                    start_point=(pos[0], line.start_point[1] + line.line_width/2 * -0.75),
                    end_point=(pos[0], line.start_point[1] + line.line_width/2 * 0.75),
                    line_width=12,
                    colour=(255, 0, 0))
                marker.plot(canvas)
                canvas.present(clear=False)

            wait_for_next(canvas)
            canvas.present(clear = False)

            # Eintrag für einen trial
            row = [exp.subject, i, "test_2", id_left, id_right, exemplars_test_1_temp[i].type, pos[0],
                   contemplation_time, rt, time, reaction_time,
                   pal_set_var, pos_set_var, bool_timepressure_first]
            writer.writerow(row)

    ####################################################################################################### Demographics
    canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
    post_test_instructions.plot(canvas)
    wait_for_next(canvas, instructions = True)

    webbrowser.open(''.join(('https://baselpsychology.eu.qualtrics.com/jfe/form/SV_0CJcBuBN2rX61aB?id=', str(exp.subject))))

expyriment.control.end()

# Balken: ein bisschen Spielraum, daneben klicken ok --> check
# falls Weiter Knopf bei Contemplation nicht rechzeitig geklickt, dann direkt nächster Trial --> check
# Familiarisierungsphsase bei Zeitdruck länger machen --> check
# in no time pressure  3 Sekunden Reaktionszeit --> check
# Familiarisierungsstimuli systematisch wählen: 111-555, 511-154,; alle Ähnlichkeitsstufen ungefähr gleich häufig --> check
# Fixation Cross kürzer --> check
# Kreis bessser platzieren
# csv.files anpassen --> check ausser counter für trials
# demographischer Fragebogen!