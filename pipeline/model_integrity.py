"""
Model artifact integrity utilities
===================================
SHA-256 verification helpers for pickled model artifacts (ASO1-2 fix).

pickle.load() can execute arbitrary code during deserialization, so
committed/distributed artifacts should be verified against a known
SHA-256 hash BEFORE loading. Expected hashes are supplied via
environment variables (see README):

  - ASO_BEST_MODEL_SHA256     -> pipeline/outputs/models/best_model.pkl
  - ASO_PREPROCESSING_SHA256  -> pipeline/outputs/models/preprocessing_pipeline.pkl

Backward compatible: if no expected hash is configured, a loud
warnings.warn is emitted and the file is loaded unverified.
"""

import hashlib
import pickle
import warnings


def sha256_of(path):
    """Compute the SHA-256 hex digest of a file's bytes (streaming)."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_pickle_verified(path, expected_sha256=None):
    """Load a pickle file, verifying its SHA-256 hash first when provided.

    Args:
        path: path to the pickle artifact.
        expected_sha256: expected SHA-256 hex digest (typically read from an
            ASO_*_SHA256 environment variable by the caller). If None, a loud
            warning is emitted and the artifact is loaded unverified
            (backward-compatible behavior).

    Raises:
        ValueError: if expected_sha256 is set and the file's actual SHA-256
            does not match (refuses to deserialize a tampered artifact).
    """
    if expected_sha256:
        actual = sha256_of(path)
        if actual.lower() != expected_sha256.strip().lower():
            raise ValueError(
                f"SHA-256 mismatch for {path}: expected {expected_sha256}, "
                f"got {actual}. Refusing to deserialize a possibly tampered artifact."
            )
    else:
        warnings.warn(
            f"Loading {path} WITHOUT integrity verification "
            "(no expected SHA-256 configured). Set the corresponding "
            "ASO_*_SHA256 environment variable to enable verification.",
            stacklevel=2,
        )
    with open(path, "rb") as f:
        return pickle.load(f)
