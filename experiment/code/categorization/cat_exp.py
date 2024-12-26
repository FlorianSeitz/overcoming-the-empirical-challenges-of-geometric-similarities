# -*- coding: utf-8 -*-
# Categorization Experiment
# 1. Import necessary libraries
import expyriment
import csv
import os
from typing import List
import numpy
import turtle
import random
import itertools
from stimulus import AdjustableStimulus
from win32api import GetSystemMetrics # Windows: install pywin32
import webbrowser

# Screen Resolution
x = GetSystemMetrics(0)
y = GetSystemMetrics(1)

# 2. Create Exemplars
# 2.1. Generate new class "Exemplar" consisting of name, values, cirt, and pos
class Exemplar:
    def __init__(self, name, values, cat, crit, pos=-1):
        self.name = name
        self.values = values
        self.cat = cat
        self.crit = crit
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
            cat = r[4]
            if i>0:
                values = (a[cue_pos_vec[0]], a[cue_pos_vec[1]], a[cue_pos_vec[2]])
                cat = r[4]
                #crit = r[5]
                ex = Exemplar(name, values, cat, -1)
                exemplars.append(ex)
            i = i+1
    return exemplars

# 2.3. Read exemplars for experiment
exemplars_training = read_csv(os.path.join(os.path.dirname(__file__), 'exemplars_training.csv'), 3)

exemplars_test = read_csv(os.path.join(os.path.dirname(__file__), 'exemplars_test.csv'), 3)
fam_blocks = 4
test_blocks = 14

left_right = expyriment.stimuli.Picture("left_right_keys.png", position = (0, 100))

# 2.4. Randomization of colour and shape
sets = list(itertools.permutations(range(3)))
pal_nr = random.sample(range(len(sets)), 1)[0]
pos_nr = random.sample(range(len(sets)), 1)[0]
# order_nr = random.sample(range(len(sets)), 1)[0]
pal_set = sets[pal_nr]
pos_set = (0, 1, 2)
pal_set_var = ''.join(str(e) for e in pal_set)

# order = sets[order_nr]
cat_order = random.sample(range(2), 2)
pos_order = random.sample(list(range(3)), 3)
pos_order_var = ''.join(str(e) for e in pos_order)

print(pal_set_var, pos_order_var)
print(cat_order)

# 3. Experiment Structure
# 3.1. Experiment initialization
exp = expyriment.design.Experiment(name="cat_exp")
expyriment.io.defaults.outputfile_time_stamp = True # omits time of experiment conduction from file name in data folder
expyriment.control.defaults.window_size = (x, y)
expyriment.control.set_develop_mode(False) # set develop mode so one does not have to render subject id, cool!
expyriment.control.initialize(exp)

# 3.2. Preload all possible stimuli
stimuli = {}
stimuli_instructions = {}
for i in numpy.array(range(1, 5)):
    for j in numpy.array(range(1, 5)):
        for k in numpy.array(range(1, 5)):
            stim = AdjustableStimulus()
            stim.init_cues(i, j, k, pal_set, pos_set, True, False)
            stimuli[str(i) + str(j) + str(k)] = stim

            stim_instructions = AdjustableStimulus(instruction="yes")
            stim_instructions.init_cues(i, j, k, pal_set, pos_set, True, False)
            stimuli_instructions[str(i) + str(j) + str(k)] = stim_instructions

# 3.2. Instructions
text_font="Century Gothic"
# 3.2.1. General instructions
general_instructions_page_1 = u"Herzlich Willkommen zu dieser Studie, in der wir untersuchen möchten, wie Menschen Kategorisierungen vornehmen. \n\n" \
                              u"Bitte lesen Sie die Instruktionen aufmerksam durch. " \
                              u"Sollten die Instruktionen unklar sein oder sollte das Experiment nicht richtig funktionieren, geben Sie bitte umgehend dem Studienleiter Bescheid. \n\n\n" \
                              u"Drücken Sie bitte die Leertaste um fortzufahren."
general_instructions_page_1 = expyriment.stimuli.TextScreen("", general_instructions_page_1, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

general_instructions_page_2a = u"In dieser Studie geht es darum zu lernen, hypothetische Produkte zu zwei Marken (Marke L oder Marke R) zuzuordnen. \n" \
                              u"Die Pfeiltaste nach links steht für Marke L. Die Pfeiltaste nach rechts steht für Marke R. \n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" \
                               u"Drücken Sie bitte die Leertaste um fortzufahren"
general_instructions_page_2a = expyriment.stimuli.TextScreen("", general_instructions_page_2a, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

general_instructions_page_2b = u"Jedes Produkt besteht aus drei Zutaten. Jedes Produkt hat eine bestimmte Menge von jeder Zutat und unterscheidet sich somit von den anderen Produkten durch eine einzigartige Kombination der Zutaten. " \
                              u"Alle Produkte bestehen aus denselben drei Zutaten, d.h. die Zutaten ändern sich nicht von Produkt zu Produkt. \n\n" \
                              u"Hier sehen Sie ein Produkt. Jede Zutat ist schematisch durch einen grauen Balken dargestellt. Die Anzahl farbiger Quadrate innerhalb der Balken gibt die Menge der jeweiligen Zutat an. \n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" \
                              u"Drücken Sie bitte die Leertaste um fortzufahren."
general_instructions_page_2b = expyriment.stimuli.TextScreen("", general_instructions_page_2b, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

general_instructions_page_3 = u"Da Sie das Experiment nicht erfolgreich beenden können, ohne die Produkte mit allen möglichen Zutaten zu kennen, sollten Sie sich jetzt die Produkte und Zutaten in Ruhe anschauen. \n\n" \
                              u"Damit Sie sich mit den Produkten vertraut machen können, klicken Sie bitte auf jede Zutat 10 mal."
general_instructions_page_3 = expyriment.stimuli.TextScreen("", general_instructions_page_3, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

number_of_clicks = u"Klicken Sie noch %i mal auf diese Zutat."

general_instructions_page_3_done = u"Weiter durch Drücken der Leertaste."
general_instructions_page_3_done = expyriment.stimuli.TextScreen("", general_instructions_page_3_done, text_size=24, background_colour=(255, 255, 255), size = (x, 0.05*y),  text_colour= (128, 128, 128), position=(0, 150), text_font=text_font)

general_instructions_page_4 = u"In jedem Durchgang wird Ihnen ein zufällig ausgewähltes Produkt mit der dazugehörigen Menge von jeder Zutat gezeigt. " \
                              u"Ihre Aufgabe ist es, richtig zu erraten, welcher Marke (Marke L oder Marke R) das jeweilige Produkt angehört, indem Sie entweder auf die linke oder die rechte Pfeiltaste drücken.\n\n\n" \
                              u"Drücken Sie bitte die Leertaste um fortzufahren."
general_instructions_page_4 = expyriment.stimuli.TextScreen("", general_instructions_page_4, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

# 3.2.1.1. Time Pressure Conditions
general_instructions_page_5_tp = u"Das Experiment gliedert sich in zwei Phasen: \n\n" \
                                 u"In der ersten Phase sehen Sie einige Produkte mehrere Male und lernen sie der richtigen Marke (Marke L oder Marke R) zuzuordnen. " \
                                 u"Nach jedem Durchgang erhalten Sie eine Rückmeldung, ob sie richtig getippt haben. " \
                                 u"Ihre Aufgabe ist zu lernen, zu welcher Marke die einzelnen Produkte gehören, und die Produkte konsistent der richtigen Marke (R oder L) zuzuordnen. " \
                                 u"Sobald Sie dies geschafft haben, gelangen Sie zur zweiten Phase. \n\n" \
                                 u"In der zweiten Phase sehen Sie erneut einige Produkte. Ihre Aufgabe ist es jedes Produkt der Marke zuzuordnen, zu der es am wahrscheinlichsten gehört (Marke L oder Marke R). " \
                                 u"In dieser Phase gibt es allerdings keine richtige Antwort mehr - Sie erhalten deswegen auch keine Rückmeldung mehr zu Ihrer Antwort. " \
                                 u"Des Weiteren gibt es in der zweiten Phase ein Zeitlimit für Ihre Antwort, das Sie nicht überschreiten sollten. \n\n\n" \
                                 u"Drücken Sie bitte die Leertaste um fortzufahren."
general_instructions_page_5_tp = expyriment.stimuli.TextScreen("", general_instructions_page_5_tp, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

general_instructions_page_5_ntp = u"Das Experiment gliedert sich in zwei Phasen: \n\n" \
                                 u"In der ersten Phase sehen Sie einige Produkte mehrere Male und lernen sie der richtigen Marke (Marke L oder Marke R) zuzuordnen. " \
                                 u"Nach jedem Durchgang erhalten Sie eine Rückmeldung, ob sie richtig getippt haben. " \
                                 u"Ihre Aufgabe ist zu lernen, zu welcher Marke die einzelnen Produkte gehören, und die Produkte konsistent der richtigen Marke (R oder L) zuzuordnen. " \
                                 u"Sobald Sie dies geschafft haben, gelangen Sie zur zweiten Phase. \n\n" \
                                 u"In der zweiten Phase sehen Sie erneut einige Produkte. Ihre Aufgabe ist es jedes Produkt der Marke zuzuordnen, zu der es am wahrscheinlichsten gehört (Marke L oder Marke R). " \
                                 u"In dieser Phase gibt es allerdings keine richtige Antwort mehr - Sie erhalten deswegen auch keine Rückmeldung mehr zu Ihrer Antwort. \n\n\n" \
                                 u"Drücken Sie bitte die Leertaste um fortzufahren."
general_instructions_page_5_ntp = expyriment.stimuli.TextScreen("", general_instructions_page_5_ntp, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

# 3.2.2. Learning Phase
train_instructions = u"Sie können nun mit der ersten Phase beginnen. \n\n" \
                     u"Drücken Sie den Pfeil nach links um ein Produkt der Marke L zuzuordnen und den Pfeil nach rechts um ein Produkt der Marke R zuzuordnen. \n\n" \
                     u"Haben Sie richtig getippt, erscheint ein freundliches Gesicht. Haben Sie falsch getippt, erscheint ein trauriges Gesicht.\n Sie können Sich weiterhin das Produkt und die Rückmeldung solange Sie möchten ansehen. " \
                     u"Durch Drücken der Pfeiltaste nach oben gelangen Sie zum nächsten Durchgang. \n\n" \
                     u"Es dauert in der Regel mehrere hundert Durchgänge, bis man anfängt, etwas über die verschiedenen Produkte zu lernen - haben Sie also bitte etwas Geduld. \n\n\n" \
                     u"Drücken Sie bitte die Leertaste um mit der ersten Phase zu beginnen."

train_instructions = expyriment.stimuli.TextScreen("", train_instructions, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

# 3.2.3. Test Phase
test_instructions_tp = u"Nun beginnt die zweite Phase des Experiments. \n\n" \
                       u"Da es in dieser Phase des Experiments keine richtigen Antworten mehr gibt, erhalten Sie nun keine Rückmeldung mehr. " \
                       u"Weisen Sie bitte jedes Produkt derjenigen Marke zu, der das Produkt für Sie am wahrscheinlichsten angehört. \n\n" \
                       u"In dieser Phase haben Sie nun zusätzlich ein Zeitlimit für jeden Durchgang. \n" \
                       u"Sie haben zwischen %d und %d Sekunden pro Durchgang Zeit. \n\n" \
                       u"Versuchen Sie bitte, dieses Zeitlimit nicht zu überschreiten! \n\n\n" \
                       u"Drücken Sie bitte die Leertaste um fortzufahren."

test_instructions_ntp = u"Nun beginnt die zweite Phase des Experiments. \n\n" \
                       u"Da es in dieser Phase des Experiments keine richtigen Antworten mehr gibt, erhalten Sie nun keine Rückmeldung mehr. " \
                       u"Weisen Sie bitte jedes Produkt derjenigen Marke zu, der das Produkt für Sie am wahrscheinlichsten angehört. \n\n\n" \
                       u"Drücken Sie bitte die Leertaste um fortzufahren."
test_instructions_ntp = expyriment.stimuli.TextScreen("", test_instructions_ntp, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

# 3.2.3.1. Familiarization
test_fam = u"Zuerst bearbeiten Sie einige Übungsdurchgänge, in denen Sie Ihnen schon bekannte Produkte sehen. \n\n" \
           u"Drücken Sie weiterhin den Pfeil nach links um ein Produkt der Marke L zuzuordnen und den Pfeil nach rechts um ein Produkt der Marke R zuzuordnen. \n\n" \
           u"Drücken Sie bitte die Leertaste um zu beginnen."
test_fam = expyriment.stimuli.TextScreen("", test_fam, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

test_test = u"Sie haben die Übungsdurchgänge abgeschlossen und werden nun weitere Produkte sehen. \n\n" \
            u"Drücken Sie weiterhin den Pfeil nach links um ein Produkt der Marke L zuzuordnen und den Pfeil nach rechts um ein Produkt der Marke R zuzuordnen. \n\n" \
            u"Drücken Sie bitte die Leertaste um zu beginnen."
test_test = expyriment.stimuli.TextScreen("", test_test, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

# 3.2.4. Post-test instruction
post_test_instructions = u"Sie haben die zweite Phase nun abgeschlossen. \n\n" \
                         u"Durch Drücken der Leertaste gelangen Sie zum demographischen Fragebogen. \n Bitte füllen Sie diesen zum Schluss noch aus."
post_test_instructions = expyriment.stimuli.TextScreen("", post_test_instructions, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, y), text_font=text_font)

# 3.2.5. Go further in instructions (spacebar)
spacebar = u"Drücken Sie die Leertaste um forzufahren."
spacebar= expyriment.stimuli.TextScreen("", spacebar, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.2*x, 0.8*y), position=(0, -0.9*y), text_font=text_font)

# 3.2.6. Feedback instructions
wrong_key = u"Verwenden Sie die Pfeiltasten: \n\n" \
            u"Drücken sie nach links für die Marke L, nach rechts für die Marke R, und nach oben um zum nächsten Durchgang zu gelangen. \n\n" \
            u"Lernen Sie Weiter durch Drücken der Pfeiltaste nach oben."
wrong_key = expyriment.stimuli.TextScreen("", wrong_key, text_bold=True, background_colour=(255, 255, 255), text_colour=(0,0,0), text_size=20, position=(0,0), size = (0.75*x, y), text_font=text_font)

too_slow_instructions = u"Leider zu langsam."
too_slow_instructions = expyriment.stimuli.TextScreen("", too_slow_instructions, text_bold=True, background_colour=(255, 255, 255), text_colour=(0,0,0), text_size=20, position=(0,-100), size = (0.75*x, y), text_font=text_font)

# 3.2.7. Feedback in training phase
incorrect_instructions = u"Falsch!"
incorrect_instructions = expyriment.stimuli.TextScreen("", incorrect_instructions, text_bold=True, background_colour=(255, 255, 255), text_colour=(255,0,0), text_size=20, position=(0,-180), size = (0.75*x, 0.1*y), text_font=text_font)

correct_instructions = u"Richtig!"
correct_instructions = expyriment.stimuli.TextScreen("", correct_instructions, text_bold=True, background_colour=(255, 255, 255), text_colour=(0,255,0), text_size=20, position=(0,-180), size = (0.75*x, 0.1*y), text_font=text_font)

# 3.2.8. Go to next trial in test phase
go_on = u"Mit dem Pfeil nach oben geht es weiter."
go_on = expyriment.stimuli.TextScreen("", go_on, text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.75*x, 0.1*y), position=(0, -300), text_font=text_font)

brand_L = u"Marke L"
brand_L  = expyriment.stimuli.TextScreen("", brand_L , text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.05*x, 0.05*y), position=(-102, -50), text_font=text_font)

brand_R = u"Marke R"
brand_R  = expyriment.stimuli.TextScreen("", brand_R , text_size=24, background_colour=(255, 255, 255), text_colour=(0, 0, 0), size = (0.05*x, 0.05*y), position=(112, -50), text_font=text_font)

# 3.3. Time Pressure Conditions
bool_timepressure = False

time_vec = []

# 3.4. Correct counter in training phase
correct = []

# 3.5. Lines for dimensions (instructions
#line1 = expyriment.stimuli.Line(start_point=(-26.25, -35), end_point=(-147.5, 35), line_width=50)
#line2 = expyriment.stimuli.Line(start_point=(26.25, -35), end_point=(147.5, 35), line_width=50)
#line3 = expyriment.stimuli.Line(start_point=(0, -80), end_point=(0, -220), line_width=50)

line1 = expyriment.stimuli.Line(start_point=(-35, -57), end_point=(-157, 13), line_width=60)
line1.rotate(15)
line2 = expyriment.stimuli.Line(start_point=(20, -11), end_point=(141, 59), line_width=60)
line2.rotate(15)
line3 = expyriment.stimuli.Line(start_point=(31, -77), end_point=(31, -217), line_width=60)
line3.rotate(15)

# 4. Blank screen
blank_screen = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)

# 5. Define writer
writer = None

# 6. Familiarization phase: Values vector
values = range(1, 5)

# 7. Fixation cross
fixation_cross = expyriment.stimuli.Shape(position = (0, 0), colour = (0, 0, 0), line_width=2)
fixation_cross.add_vertices([(25, 0), (-50, 0), (25, 0), (0, 25), (0, -50), (0, 25)])

time_pressure_in_ms = "NA"
####################################################################################################### Start Experiment
expyriment.control.start()
with open(str("data/results_" + str(exp.subject) + ".csv"), "wb") as csvfile:
    writer = csv.writer(csvfile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    # Header
    row = ["subj_id", "trial", "block", "stim", "true_cat", "cat0_is_right_key", "color_presentation_order", "feature_presentation_order", "response", "time", "time_pressure_cond", "time_pressure_in_ms", "too_slow"]
    writer.writerow(row)

    if (exp.subject % 2) == 1:
        bool_timepressure = True

    canvas = expyriment.stimuli.Canvas(size=(1920, 1080), colour=expyriment.misc.constants.C_WHITE)

    ############################################################################################## Familiarization Phase
    # Instructions
    general_instructions_page_1.plot(canvas)
    canvas.present()
    exp.clock.wait(3000)
    exp.keyboard.wait(expyriment.misc.constants.K_SPACE)

    general_instructions_page_2a.plot(canvas)
    left_right.plot(canvas)
    brand_L.plot(canvas)
    brand_R.plot(canvas)
    canvas.present()
    exp.clock.wait(3000)
    exp.keyboard.wait(expyriment.misc.constants.K_SPACE)

    start_vec = (random.sample(range(1, 5), 3))
    stim = stimuli_instructions[''.join(str(e) for e in start_vec)]

    general_instructions_page_2b.plot(canvas)
    stim.canvas_stim.plot(canvas)
    canvas.present()
    exp.clock.wait(3000)
    exp.keyboard.wait(expyriment.misc.constants.K_SPACE)

    # Dimension experiencing
    exp.mouse.show_cursor()
    clicks_necessary = list((10, 10, 10))
    direction = list((1, 1, 1))
    for i in list((0, 1, 2)):
        if start_vec[i] == 4:
            direction[i] = -1

    while True:
        canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
        general_instructions_page_3.plot(canvas)
        if numpy.sum(clicks_necessary) == 0: # and keyboard.is_pressed('a'):
            general_instructions_page_3_done.plot(canvas)
        stim = stimuli_instructions[''.join(str(e) for e in start_vec)]
        stim.canvas_stim.plot(canvas)

        number_of_clicks1 = number_of_clicks % clicks_necessary[0]
        number_of_clicks2 = number_of_clicks % clicks_necessary[1]
        number_of_clicks3 = number_of_clicks % clicks_necessary[2]
        number_of_clicks1 = expyriment.stimuli.TextScreen("", number_of_clicks1, text_size=18,
                                                         background_colour=(255, 255, 255), text_colour=(0, 0, 0),
                                                         text_bold=True, size=(0.2 * x, 0.03*y), position=(-345, 50), text_font=text_font)
        number_of_clicks2 = expyriment.stimuli.TextScreen("", number_of_clicks2, text_size=18,
                                                         background_colour=(255, 255, 255), text_colour=(0, 0, 0),
                                                         text_bold=True, size=(0.2 * x, 0.03*y), position = (345, 50), text_font=text_font)
        number_of_clicks3 = expyriment.stimuli.TextScreen("", number_of_clicks3, text_size=18,
                                                         background_colour=(255, 255, 255), text_colour=(0, 0, 0),
                                                         text_bold=True, size=(0.2 * x, 0.03*y), position = (0, -260), text_font=text_font)
        number_of_clicks1.plot(canvas)
        number_of_clicks2.plot(canvas)
        number_of_clicks3.plot(canvas)

        canvas.present()
        if numpy.sum(clicks_necessary) == 0:
            break

        click = exp.mouse.wait_press(0)
        pos = click[1]

        index = 999

        if line1.overlapping_with_position(pos): index = 0
        if line2.overlapping_with_position(pos): index = 1
        if line3.overlapping_with_position(pos): index = 2
        if index != 999:
            start_vec[index] = start_vec[index] + direction[index]
            if clicks_necessary[index] > 0:
                clicks_necessary[index] = clicks_necessary[index] - 1
        for i in list((0, 1, 2)):
            if start_vec[i] == 4:
                direction[i] = -1
            if start_vec[i] == 1:
                direction[i] = 1

    exp.keyboard.wait(expyriment.misc.constants.K_SPACE)
    exp.mouse.hide_cursor()

    # Instructions Part 2
    general_instructions_page_4.plot(canvas)
    canvas.present()
    exp.clock.wait(3000)
    exp.keyboard.wait(expyriment.misc.constants.K_SPACE)

    if bool_timepressure:
        general_instructions_page_5_tp.plot(canvas)
    else:
        general_instructions_page_5_ntp.plot(canvas)
    canvas.present()
    exp.clock.wait(3000)
    exp.keyboard.wait(expyriment.misc.constants.K_SPACE)

    ##################################################################################################### Training Phase
    # Instructions
    train_instructions.plot(canvas)
    canvas.present()
    exp.clock.wait(3000)
    exp.keyboard.wait(expyriment.misc.constants.K_SPACE)

    # Stimulus presentation
    trial_counter = 0
    accuracy_achieved = False
    while accuracy_achieved == False:
        if trial_counter % 8 == 0:
            random_order_training = numpy.array(range(len(exemplars_training)))
            random_order_training = random.sample(random_order_training, len(random_order_training))

        for i in random_order_training:
            blank_screen.present()
            exp.clock.wait(1000)

            key = 0
            canvas = expyriment.stimuli.Canvas(size=(x, y), colour=expyriment.misc.constants.C_WHITE)
            name = exemplars_training[i].name
            name = ''.join(str(e) for e in numpy.array([name[pos_order[0]], name[pos_order[1]], name[pos_order[2]]]))
            stim = stimuli[name]
            stim.canvas_stim.plot(canvas)
            while key != 275 and key != 276:
                canvas.present()
                key, rt = exp.keyboard.wait()
                if key != 275 and key != 276:
                    wrong_key.plot(canvas)
                    canvas.present()
                    exp.keyboard.wait(expyriment.misc.constants.K_UP)
                    canvas = expyriment.stimuli.Canvas(size=(x, y), colour=expyriment.misc.constants.C_WHITE)
                    stim.canvas_stim.plot(canvas)
            if key == 276: response = str(cat_order[0]) # 276 = left
            if key == 275: response = str(cat_order[1])

            if trial_counter < 2:
                go_on.plot(canvas)
            if response == exemplars_training[i].cat:
                correct_instructions.plot(canvas)
                canvas.present(clear=False, update=False)
                head = expyriment.stimuli.Circle(25, position=(0, -100), colour=(0,255,0))
                head.present(clear=False, update = False)
                mouth = expyriment.stimuli.Circle(16, position=(0, -104), colour=(0,0,0))
                mouth.present(clear=False, update = False)
                mouth = expyriment.stimuli.Rectangle(size=(32, 16), colour=(0,255,0), position= (0, -96))
                mouth.present(clear=False, update = False)
                eye = expyriment.stimuli.Circle(2, position=(-5, -88), colour=(0,0,0))
                eye.present(clear=False, update = False)
                eye = expyriment.stimuli.Circle(2, position=(5, -88), colour=(0,0,0))
                eye.present(clear=False)

                correct.append(1)

            else:
                incorrect_instructions.plot(canvas)
                canvas.present(clear=False, update=False)
                head = expyriment.stimuli.Circle(25, position=(0, -100), colour=(255,0,0))
                head.present(clear=False, update = False)
                mouth = expyriment.stimuli.Ellipse((16,8), position=(0, -112), colour=(0,0,0))
                mouth.present(clear=False, update = False)
                mouth = expyriment.stimuli.Rectangle(size=(32, 8), colour=(255,0,0), position= (0, -116))
                mouth.present(clear=False, update = False)
                eye = expyriment.stimuli.Circle(2, position=(-5, -88), colour=(0,0,0))
                eye.present(clear=False, update = False)
                eye = expyriment.stimuli.Circle(2, position=(5, -88), colour=(0,0,0))
                eye.present(clear=False)

                correct.append(0)

            time_vec.append(rt)
            trial_counter = trial_counter + 1

            # Eintrag für einen trial
            row = [exp.subject, trial_counter, "training", exemplars_training[i].name, exemplars_training[i].cat, cat_order[0], pal_set_var, pos_order_var, response, rt, bool_timepressure, time_pressure_in_ms, False]
            writer.writerow(row)

            if numpy.mean(correct[-100:]) >= .80 and numpy.mean(correct[-24:]) == 1 and trial_counter > 100:
                accuracy_achieved = True
                break

            if (trial_counter % 50) == 0 and trial_counter != 50:
                exp.keyboard.wait(expyriment.misc.constants.K_UP)

                accuracy = u"Momentan klassifizieren Sie %d Prozent der Produkte richtig. Weiter so! \n\n" \
                           u"Versuchen Sie weiterzulernen, um noch mehr richtig zu klassifizieren. \n\n" \
                           u"Mit der Pfeiltaste nach oben geht es mit dem Lernen weiter." % (numpy.mean(correct[-100:]) * 100.0)
                accuracy = expyriment.stimuli.TextScreen("", accuracy, text_size=24, background_colour=(255, 255, 255),
                                                       text_colour=(0, 0, 0), size=(0.75 * x, y), text_font=text_font)
                accuracy.plot(canvas)
                canvas.present()

            exp.keyboard.wait(expyriment.misc.constants.K_UP)

    # Determine Time Pressure
    if bool_timepressure:
		time = 400 + (numpy.median(time_vec[-100:]) * .3)
		time_pressure_in_ms = numpy.round(time)

    ######################################################################################################### Test Phase
    # Instructions
    if bool_timepressure:
        test_instructions_tp = test_instructions_tp % (numpy.floor(time / 1000.0), numpy.ceil(time / 1000.0))
        test_instructions_tp = expyriment.stimuli.TextScreen("", test_instructions_tp, text_size=24,
                                                             background_colour=(255, 255, 255), text_colour=(0, 0, 0),
                                                             size=(0.75 * x, y), text_font=text_font)
        test_instructions_tp.plot(canvas)
    else:
        test_instructions_ntp.plot(canvas)
    canvas.present()
    exp.clock.wait(3000)
    exp.keyboard.wait(expyriment.misc.constants.K_SPACE)

    test_fam.plot(canvas)
    canvas.present()
    exp.clock.wait(3000)
    exp.keyboard.wait(expyriment.misc.constants.K_SPACE)

    # Stimulus Presentation
    blank_screen.present()
    stim = AdjustableStimulus()

    # Habituation
    for block in range(0, fam_blocks):
        random_order_training = numpy.array(range(len(exemplars_training)))
        random_order_training = random.sample(random_order_training, len(random_order_training))
 
        for i in random_order_training:
            canvas = expyriment.stimuli.Canvas(size=(x, y), colour=expyriment.misc.constants.C_WHITE)

            name = exemplars_training[i].name
            name = ''.join(str(e) for e in numpy.array([name[pos_order[0]], name[pos_order[1]], name[pos_order[2]]]))
            stim = stimuli[name]
            stim.canvas_stim.plot(canvas)
            canvas.present()
            if bool_timepressure:
                key, rt = exp.keyboard.wait({expyriment.misc.constants.K_LEFT, expyriment.misc.constants.K_RIGHT},
                                            duration=time)
                canvas = expyriment.stimuli.Canvas(size=(x, y), colour=expyriment.misc.constants.C_WHITE)
                if key is None:
                    too_slow = True
                    too_slow_instructions.plot(canvas)
                    canvas.present()
                    key, rt = exp.keyboard.wait({expyriment.misc.constants.K_LEFT, expyriment.misc.constants.K_RIGHT},
                                                duration=1000)
                    print rt
                    if key is None:
                        rt = "NA"
                        key = "NA"
                        response = "NA"
                    else:
                        exp.clock.wait(1000-rt)
                        rt = time + rt
                else:
                    too_slow = False

            else:
                too_slow = False
                key, rt = exp.keyboard.wait({expyriment.misc.constants.K_LEFT, expyriment.misc.constants.K_RIGHT})
                canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
            if key == 276: response = str(cat_order[0]) # 276 = left
            if key == 275: response = str(cat_order[1])
            trial_counter = trial_counter + 1

            row = [exp.subject, trial_counter, "familiarization", exemplars_training[i].name, exemplars_training[i].cat, cat_order[0], pal_set_var, pos_order_var, response, rt, bool_timepressure, time_pressure_in_ms, too_slow]
            writer.writerow(row)

            canvas = expyriment.stimuli.Canvas(size=(x, y), colour=expyriment.misc.constants.C_WHITE)
            fixation_cross.plot(canvas)
            canvas.present()
            exp.clock.wait(500)

    # Critical blocks
    test_test.plot(canvas)
    canvas.present()
    exp.clock.wait(3000)
    exp.keyboard.wait(expyriment.misc.constants.K_SPACE)

    for block in range(0, test_blocks):
        random_order_test = numpy.array(range(len(exemplars_test)))
        random_order_test = random.sample(random_order_test, len(random_order_test))

        for i in random_order_test:
            canvas = expyriment.stimuli.Canvas(size=(x, y), colour=expyriment.misc.constants.C_WHITE)

            name = exemplars_test[i].name
            name = ''.join(str(e) for e in numpy.array([name[pos_order[0]], name[pos_order[1]], name[pos_order[2]]]))
            stim = stimuli[name]
            stim.canvas_stim.plot(canvas)
            canvas.present()

            if bool_timepressure:
                key, rt = exp.keyboard.wait({expyriment.misc.constants.K_LEFT, expyriment.misc.constants.K_RIGHT},
                                            duration=time)
                canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
                if key is None:
                    too_slow = True
                    too_slow_instructions.plot(canvas)
                    canvas.present()
                    key, rt = exp.keyboard.wait({expyriment.misc.constants.K_LEFT, expyriment.misc.constants.K_RIGHT},
                                                duration=1000)
                    if key is None:
                        rt = "NA"
                        key = "NA"
                        response = "NA"
                    else:
                        exp.clock.wait(1000-rt)
                        rt = time + rt
                else:
                    too_slow = False

            else:
                too_slow = False
                key, rt = exp.keyboard.wait({expyriment.misc.constants.K_LEFT, expyriment.misc.constants.K_RIGHT})
                canvas = expyriment.stimuli.BlankScreen(colour=expyriment.misc.constants.C_WHITE)
            if key == 276: response = str(cat_order[0]) # 276 = left
            if key == 275: response = str(cat_order[1])
            trial_counter = trial_counter + 1
			
			# Eintrag für einen trial
            row = [exp.subject, trial_counter, "test", exemplars_test[i].name, "NA", cat_order[0], pal_set_var, pos_order_var, response, rt, bool_timepressure, time_pressure_in_ms, too_slow]
            writer.writerow(row)

            canvas = expyriment.stimuli.Canvas(size=(x, y), colour=expyriment.misc.constants.C_WHITE)
            fixation_cross.plot(canvas)
            canvas.present()
            exp.clock.wait(500)

    ####################################################################################################### Demographics
    post_test_instructions.plot(canvas)
    canvas.present()
    exp.clock.wait(3000)
    exp.keyboard.wait(expyriment.misc.constants.K_SPACE)

    webbrowser.open(''.join(('https://baselpsychology.eu.qualtrics.com/jfe/form/SV_3f5jiRO7CLOHc2x?id=', str(exp.subject))))

expyriment.control.end()

# Notes:
# Den Probanden klar machen, dass sie beide Kategoriestrukturen lernen müssen, nicht nur eine
# Testphase "Mit der Pfeiltaste nach oben geht es weiter."
# Testphase: sagen dass zwei Familiarisierungsblöcke
# unschön: Ende der Familiarisierung, man kann nicht mehr weiterklicken