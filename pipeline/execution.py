"""Deterministic execution helpers that do not alter scientific parameters."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
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


def ordered_process_map(
    function: Callable[[Input], Output],
    items: Iterable[Input],
    *,
    workers: int,
    initializer: Callable[..., None] | None = None,
    initargs: tuple[object, ...] = (),
) -> list[Output]:
    """Map CPU-bound jobs in input order using initialized worker processes."""

    materialized = list(items)
    if workers < 1:
        raise ValueError("workers must be at least one")
    if workers == 1 or len(materialized) < 2:
        if initializer is not None:
            initializer(*initargs)
        return [function(item) for item in materialized]
    with ProcessPoolExecutor(
        max_workers=workers, initializer=initializer, initargs=initargs
    ) as executor:
        return list(executor.map(function, materialized))
