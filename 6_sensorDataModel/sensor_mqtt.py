import paho.mqtt.client as mqtt
import pandas as pd
import time
import os

print(os.getcwd())
def sget_sensor_mqtt():
    test_data = {'RHumidity': [70.00],
    'RTemperature (C)': [26.70],
    'RTemperature (F)': [80.06],
    'WHumidity': [16.8],
    'WTemperature': [91],
    'WPressure': [1005],
    'Timestamp': ["2023-07-31 14:01:23,62.00"]
    }
    df = pd.DataFrame(test_data)
    df.to_csv("ye1.csv", index=False)
    return df


def get_sensor_mqtt():
    
    # Counter to track number of received messages
    message_count = 0

    # Callback function triggered when client connects to the broker
    def on_connect(client, userdata, flags, rc):
        print("Connected to broker")
        # Subscribe to the desired MQTT topic
        client.subscribe("sensor_weather_data")

    # Callback function triggered when a new message is received
    def on_message(client, userdata, msg):
        nonlocal message_count
        global df
        data = {"Timestamp":[], "RHumidity":[], "RTemperature (C)":[], "RTemperature (F)":[], "WHumidity":[], "WTemperature":[], "WPressure":[], "Weather Conditions":[]}
        df = pd.DataFrame(data)
        print(df.shape)
        print(f"Received message: {msg.payload.decode()}")
        # Split the payload into separate values
        values = msg.payload.decode().split(',')

        # Create a new DataFrame with the values
        new_df = pd.DataFrame([values], columns=data)

        # Concatenate the new DataFrame with the existing df DataFrame
        df = pd.concat([df, new_df], ignore_index=True)
  
        # Increment the message count
        message_count += 1
        
        # Stop the MQTT loop after receiving 10 messages
        if message_count >= 1:
            client.disconnect()        
        
        return df

    # Create a MQTT client instance
    client = mqtt.Client()

    # Assign the callback functions
    client.on_connect = on_connect
    client.on_message = on_message

    # Connect to the MQTT broker
    client.connect("mqtt.eclipseprojects.io", 1883, 60)

    # Start the MQTT loop in a blocking mode
    client.loop_forever()

    return df

