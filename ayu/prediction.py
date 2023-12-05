import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import argparse
import pickle
import pandas as pd
import os
import xgboost as xgb

import ayu.preprocessing

try:
    from importlib import resources as impresources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as impresources

from . import models


def load_xgboost_model(model_file):
    model = xgb.XGBClassifier()
    model.load_model(model_file)
    return model

def load_predictor_models():
    predictor_dict = {
        'cyto':(impresources.files(models) / 'xgb_cyto.json'),
        'peri':(impresources.files(models) / 'xgb_peri.json'),
        'extr':(impresources.files(models) / 'xgb_extr.json'),
    }

    for key in predictor_dict:
        predictor_dict[key] = load_xgboost_model(predictor_dict[key])
    
    return predictor_dict


def run_predictions(feature_file, predictor_dict):
    to_predict_df_raw = pd.read_csv(feature_file, sep='\t', comment='#')

    finished_predictions_dict = {}
    for key in predictor_dict:
        iteration_end = int(predictor_dict[key].get_booster().attributes()['best_iteration'])
        cols_when_model_builds = predictor_dict[key].get_booster().feature_names
        to_predict_df = to_predict_df_raw[cols_when_model_builds]
        finished_predictions_dict[key] = predictor_dict[key].predict_proba(to_predict_df, iteration_range=(0, iteration_end))
    
    results_dict = {'prot_ID':to_predict_df_raw['prot_ID'].values.tolist()}
    for key in finished_predictions_dict:
        results_dict[key] = [x[1] for x in finished_predictions_dict[key]]
    
    results_df = pd.DataFrame(results_dict)
    
    return results_df

def save_prediction_file(result_df, out_file):
    result_df = result_df.reindex(sorted(result_df.columns), axis=1)
    result_df.to_csv(out_file, mode='a', header=not os.path.exists(out_file), index=False, sep = '\t')

