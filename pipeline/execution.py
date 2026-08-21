"""Deterministic execution helpers that do not alter scientific parameters."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from collections.abc import Callable, Iterable
from typing import TypeVar


Input = TypeVar("Input")
Output = TypeVar("Output")


def ordered_thread_map(
    function: Callable[[Input], Output],
    items: Iterable[Input],
    *,
    workers: int,
) -> list[Output]:
    """Map independent jobs concurrently while retaining exact input order."""

    materialized = list(items)
    if workers < 1:
        raise ValueError("workers must be at least one")
    if workers == 1 or len(materialized) < 2:
        return [function(item) for item in materialized]
    with ThreadPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(function, materialized))
