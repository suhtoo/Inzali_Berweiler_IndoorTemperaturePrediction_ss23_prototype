import pandas as pd
import numpy as np
import random
import os

os.chdir('/Users/hiz/Desktop/RWTH_3rdSem/CR_Protyping Project DSS - BIM/Project/Inzali_Berweiler_IndoorTemperaturePrediction_ss23_prototype/4_deploymentModule')

# Function to generate random data for the DataFrame
def generate_random_data():
    windspeed = random.uniform(0, 10)
    sunHeight = random.uniform(0, 90)
    weatherTemp = random.uniform(273, 320)
    windowStatus = random.choice([1, 2, 3])
    return windspeed, sunHeight, weatherTemp, windowStatus

# Function to create a DataFrame with random data and datetime index
def create_random_dataframe(start_time, end_time):
    index = pd.date_range(start=start_time, end=end_time, freq='10S')
    data = [generate_random_data() for _ in range(len(index))]
    df = pd.DataFrame(data, columns=['windspeed', 'sunHeight', 'weatherTemp', 'windowStatus'], index=index)
    return df

if __name__ == "__main__":
    start_time = pd.Timestamp.now()
    end_time = start_time + pd.Timedelta(minutes=10)
    random_df = create_random_dataframe(start_time, end_time)

    # Save the DataFrame to a CSV file
    random_df.to_csv("random_data.csv", index=True)
    print("random_data.csv'")
