import PySimpleGUI as sg
import predictions
import time
layout = [[sg.Text('steps:')],
          [sg.Input(key='-INPUT-', default_text="100", size=(30, 5))],
          [sg.Text('predicted Temperature:')],
          [sg.Text('', key='-OUTPUT-')],
          [sg.Button('Start')], [sg.Button('Stop')]]

window = sg.Window('My Window', layout)

while True:
    event, values = window.read()

    time.sleep(1)
    # Initialize flag variable
    continue_program = True
    
    if event == sg.WINDOW_CLOSED:
        break
    if event == 'Start':
        n = values.get('-INPUT-', '')
        try:
            n = int(n)  # Convert the input value to an integer

            if 'current_data_test' in dir(predictions):
                while continue_program:  # Loop until continue_program is False
                    value_l = predictions.current_data_test(n)
                    if value_l:
                        value = value_l[0]
                        formatted_value = f"{value}°C"
                        window['-OUTPUT-'].update(formatted_value)
                    else:
                        sg.popup('Error: Failed to retrieve value from predictions.py')
                    
                    # Check if the Stop button is pressed
                    event, _ = window.read(timeout=100)  # Check for events every 100 milliseconds
                    if event == 'Stop':
                        continue_program = False  # Set continue_program to False

            else:
                sg.popup('Error: Function current_data_test not found in predictions.py')
        except ValueError:
            sg.popup('Error: Please enter a valid integer')

window.close()

