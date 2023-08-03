import joblib
import pandas as pd
import os

os.chdir('/Users/hiz/Desktop/RWTH_3rdSem/CR_Protyping Project DSS - BIM/Project/Inzali_Berweiler_IndoorTemperaturePrediction_ss23_prototype_local/4_deploymentModule')

# Load the trained ML model
rf_model = joblib.load("/Users/hiz/Desktop/RWTH_3rdSem/CR_Protyping Project DSS - BIM/Project/Inzali_Berweiler_IndoorTemperaturePrediction_ss23_prototype_local/3_modelTrainingModule/rf_model.joblib")

# Function to make roomTemp predictions using the model and the DataFrame
def make_roomTemp_predictions(model, dataframe):
    roomTemp_predictions = model.predict(dataframe)
    return roomTemp_predictions

if __name__ == "__main__":
    # Load the DataFrame from the CSV file
    random_df = pd.read_csv("random_data.csv", index_col=0)
    roomTemp_predictions = make_roomTemp_predictions(rf_model, random_df)

    # Concatenate the original DataFrame with the predictions
    predictions_df = random_df.copy()
    predictions_df["roomTemp_predictions"] = roomTemp_predictions

    # Print the predictions to the console
    print(predictions_df)

    # Save the predictions to a file (e.g., CSV)
    predictions_df.to_csv("predictions.csv")

    print("Room temperature predictions with original data saved to 'predictions.csv'")
