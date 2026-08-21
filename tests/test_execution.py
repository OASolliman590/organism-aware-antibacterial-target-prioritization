import pytest

from pipeline.execution import ordered_thread_map


def test_ordered_thread_map_preserves_order_and_propagates_errors():
    assert ordered_thread_map(lambda value: value * value, [3, 1, 2], workers=3) == [9, 1, 4]

    def fail_on_two(value):
        if value == 2:
            raise RuntimeError("declared worker failure")
        return value

    with pytest.raises(RuntimeError, match="declared worker failure"):
        ordered_thread_map(fail_on_two, [1, 2, 3], workers=2)
