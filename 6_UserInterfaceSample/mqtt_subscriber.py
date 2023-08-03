import paho.mqtt.client as mqtt
import csv

# MQTT Broker details
broker_address = "mqtt.eclipseprojects.io"  # MQTT broker address
broker_port = 1883  # MQTT broker port
topic = "room_environment_data"  # Topic to subscribe to

# Column names for the CSV file
column_names = ["timestamp", "windspeed", "sunHeight", "weatherTemp", "windowStatus", "roomTemp_predictions"]

# Callback function when the subscriber connects to the MQTT broker
def on_connect(client, userdata, flags, rc):
    if rc == 0:
        print("Connected to MQTT Broker!")
        client.subscribe(topic)
    else:
        print("Failed to connect, return code %d\n", rc)

# Callback function when a message is received from the MQTT broker
def on_message(client, userdata, msg):
    received_data = msg.payload.decode("utf-8").split(',')
    print(f"Received data: {received_data}")
    
    # Save received data to CSV file
    with open("received_data.csv", "a", newline="") as csvfile:
        writer = csv.writer(csvfile)
        # Check if the CSV file is empty to write the header row
        if csvfile.tell() == 0:
            writer.writerow(column_names)
        writer.writerow(received_data)

# Create MQTT client instance and set callback functions
client = mqtt.Client()
client.on_connect = on_connect
client.on_message = on_message

# Connect to MQTT broker and start the loop
client.connect(broker_address, broker_port)
client.loop_forever()
