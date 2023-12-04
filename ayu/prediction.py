import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import argparse
import pickle
import pandas as pd
import os
import xgboost as xgb

def load_xgboost_model(model_file):
    model = xgb.XGBClassifier()
    model.load_model(model_file)
    return model

