import pandas as pd
import numpy as np
import warnings
import pickle
warnings.filterwarnings("ignore")

from sklearn.metrics import precision_score, recall_score
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score, classification_report
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
import optuna
from optuna.samplers import TPESampler
optuna.logging.set_verbosity(optuna.logging.WARNING)

df = pd.read_csv("features.csv")
X  = df.iloc[:, 1:-1].values
y  = df.iloc[:, -1].values
feature_names = np.array(df.columns[1:-1].tolist())

print(f"Dataset: {X.shape[0]} samples, {X.shape[1]} features")
print(f"Class distribution: {dict(zip(*np.unique(y, return_counts=True)))}")

X_train_full, X_test_full, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=413, stratify=y)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=413)

imputer = SimpleImputer(strategy="median")
X_train_full = imputer.fit_transform(X_train_full)
X_test_full  = imputer.transform(X_test_full)

imp_df = pd.read_csv("feature_importance_rf.csv")
ranked_features = imp_df["feature"].tolist()

feat_to_idx = {name: idx for idx, name in enumerate(feature_names)}

summary_rows = []

for n_top in [10, 20, 30, 40, 50]:
    print(f"\n{'='*50}")
    print(f"  Top-{n_top} features")
    print(f"{'='*50}")

    top_feats   = ranked_features[:n_top]
    top_idx     = [feat_to_idx[f] for f in top_feats if f in feat_to_idx]
    X_train     = X_train_full[:, top_idx]
    X_test      = X_test_full[:,  top_idx]

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

    study = optuna.create_study(direction="maximize", sampler=TPESampler(seed=42))
    study.optimize(rf_objective, n_trials=200, show_progress_bar=True)

    print(f"  Best CV AUC : {study.best_value:.4f}")
    print(f"  Best params : {study.best_params}")

    best_rf = RandomForestClassifier(**study.best_params, random_state=413, n_jobs=-1)
    best_rf.fit(X_train, y_train)

    y_pred = best_rf.predict(X_test)
    y_prob = best_rf.predict_proba(X_test)[:, 1]

    print(f"  CV AUC  : {study.best_value:.4f}")
    print(f"  Test AUC: {roc_auc_score(y_test, y_prob):.4f}")
    print(f"  Accuracy: {accuracy_score(y_test, y_pred):.4f}")
    print(f"  F1      : {f1_score(y_test, y_pred):.4f}")
    print(classification_report(y_test, y_pred))

    pd.DataFrame({
        "feature":    top_feats,
        "importance": best_rf.feature_importances_,
    }).sort_values("importance", ascending=False).to_csv(
        f"feature_importance_rf_{n_top}.csv", index=False)

    with open(f"rf_model_{n_top}.pkl", "wb") as f:
        pickle.dump({"model": best_rf, "imputer": imputer,
                     "top_features": top_feats, "top_idx": top_idx}, f)

    train_prob = best_rf.predict_proba(X_train)[:, 1]
    train_pred = best_rf.predict(X_train)

    pd.DataFrame({"true": y_train, "pred": train_pred, "prob": train_prob}
                 ).to_csv(f"train_predictions_{n_top}.csv", index=False)

    pd.DataFrame({"true": y_test, "pred": y_pred, "prob": y_prob}
                 ).to_csv(f"test_predictions_{n_top}.csv", index=False)

    for split, X_sp, y_sp, y_pr, y_pb in [
        ("Training set", X_train, y_train, train_pred, train_prob),
        ("Test set",     X_test,  y_test,  y_pred,     y_prob),
    ]:
        summary_rows.append({
            "Dataset":   f"Top-{n_top}",
            "Split":     split,
            "Acc":       round(accuracy_score(y_sp, y_pr), 4),
            "Prec":      round(precision_score(y_sp, y_pr), 4),
            "Recall":    round(recall_score(y_sp, y_pr), 4),
            "F1":        round(f1_score(y_sp, y_pr), 4),
            "AUC":       round(roc_auc_score(y_sp, y_pb), 4),
        })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv("model_selection_summary.csv", index=False)
print(summary_df.to_string(index=False))