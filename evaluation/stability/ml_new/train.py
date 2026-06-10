import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.metrics import (accuracy_score, roc_auc_score,
                             f1_score, classification_report)
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
import optuna
from optuna.samplers import TPESampler
optuna.logging.set_verbosity(optuna.logging.WARNING)

df = pd.read_csv("features.csv")

X  = df.iloc[:, 1:-1].values
y  = df.iloc[:, -1].values
feature_names = df.columns[1:-1].tolist()

print(f"Dataset: {X.shape[0]} samples, {X.shape[1]} features")
print(f"Class distribution: {dict(zip(*np.unique(y, return_counts=True)))}")

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=413, stratify=y)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=413)

imputer = SimpleImputer(strategy="median")
X_train = imputer.fit_transform(X_train)
X_test  = imputer.transform(X_test)

def rf_objective(trial):
    params = {
        "n_estimators":      trial.suggest_int("n_estimators", 50, 1000),
        "max_depth":         trial.suggest_int("max_depth", 3, 50),
        "min_samples_split": trial.suggest_int("min_samples_split", 2, 20),
        "min_samples_leaf":  trial.suggest_int("min_samples_leaf", 1, 10),
        "max_features":      trial.suggest_categorical("max_features", ["sqrt", "log2", 0.2, 0.3, 0.5]),
        "bootstrap":         trial.suggest_categorical("bootstrap", [True, False]),
        "random_state": 413,
        "n_jobs": -1,
    }
    model = RandomForestClassifier(**params)
    return cross_val_score(model, X_train, y_train, cv=cv,
                           scoring="roc_auc", n_jobs=-1).mean()

print("\nOptimizing Random Forest...")
study = optuna.create_study(direction="maximize", sampler=TPESampler(seed=42))
study.optimize(rf_objective, n_trials=200, show_progress_bar=True)

print(f"\n  Best CV AUC : {study.best_value:.4f}")
print(f"  Best params : {study.best_params}")

best_rf = RandomForestClassifier(**study.best_params, random_state=413, n_jobs=-1)
best_rf.fit(X_train, y_train)

y_pred = best_rf.predict(X_test)
y_prob = best_rf.predict_proba(X_test)[:, 1]

print(f"\n{'='*50}")
print(f"  CV AUC  : {study.best_value:.4f}")
print(f"  Test AUC: {roc_auc_score(y_test, y_prob):.4f}")
print(f"  Accuracy: {accuracy_score(y_test, y_pred):.4f}")
print(f"  F1      : {f1_score(y_test, y_pred):.4f}")
print(f"\n{classification_report(y_test, y_pred)}")

imp_df = pd.DataFrame({
    "feature":    feature_names,
    "importance": best_rf.feature_importances_
}).sort_values("importance", ascending=False)

imp_df.to_csv("feature_importance_rf.csv", index=False)

import pickle

with open("rf_model.pkl", "wb") as f:
    pickle.dump({"model": best_rf, "imputer": imputer}, f)

train_prob = best_rf.predict_proba(X_train)[:, 1]
train_pred = best_rf.predict(X_train)

pd.DataFrame({
    "experiment":  y_train,
    "pred":  train_pred,
    "prob":  train_prob,
}).to_csv("train_predictions.csv", index=False)

pd.DataFrame({
    "experiment": y_test,
    "pred": y_pred,
    "prob": y_prob,
}).to_csv("test_predictions.csv", index=False)
