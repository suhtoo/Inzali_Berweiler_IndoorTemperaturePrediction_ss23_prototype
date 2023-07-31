import joblib
import pandas as pd
import os

os.chdir('/Users/hiz/Desktop/RWTH_3rdSem/CR_Protyping Project DSS - BIM/Project/Inzali_Berweiler_IndoorTemperaturePrediction_ss23_prototype/4_deploymentModule')

# Load the trained ML model
rf_model = joblib.load("/Users/hiz/Desktop/RWTH_3rdSem/CR_Protyping Project DSS - BIM/Project/Inzali_Berweiler_IndoorTemperaturePrediction_ss23_prototype/3_modelTrainingModule/rf_model.joblib")

# Function to make roomTemp predictions using the model and the DataFrame
def make_roomTemp_predictions(model, dataframe):
    roomTemp_predictions = model.predict(dataframe)
    return roomTemp_predictions

if __name__ == "__main__":
    # Load the DataFrame from the CSV file
    random_df = pd.read_csv("random_data.csv", index_col=0)
    roomTemp_predictions = make_roomTemp_predictions(rf_model, random_df)
    print(roomTemp_predictions)

    # Save the predictions to a file (e.g., CSV)
    pd.DataFrame({"roomTemp_predictions": roomTemp_predictions}).to_csv("predictions.csv", index=False)
    print("Room temperature predictions saved to 'predictions.csv'")