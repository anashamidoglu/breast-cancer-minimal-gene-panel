"""Reusable estimators with fold-local class balancing."""

from sklearn.utils.class_weight import compute_sample_weight
from xgboost import XGBClassifier


class BalancedXGBClassifier(XGBClassifier):
    """Compute class-balancing weights from the labels supplied to each fit."""

    def fit(self, X, y, **kwargs):
        if "sample_weight" in kwargs:
            raise ValueError("Weights must be computed from the current fit labels only.")
        return super().fit(X, y, sample_weight=compute_sample_weight("balanced", y), **kwargs)
