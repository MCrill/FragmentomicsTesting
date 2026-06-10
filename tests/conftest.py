import sys
from pathlib import Path

import pytest

# make the local tests/ helpers importable as a flat module
sys.path.insert(0, str(Path(__file__).parent))

from synthetic import make_reference, simulate_fragments  # noqa: E402


@pytest.fixture
def reference():
    return make_reference()


@pytest.fixture
def biased_sample(reference):
    frags, edited_ref, dyads = simulate_fragments(reference, phase_bias=True)
    return frags, edited_ref, dyads


@pytest.fixture
def unbiased_sample(reference):
    frags, edited_ref, dyads = simulate_fragments(reference, phase_bias=False)
    return frags, edited_ref, dyads
