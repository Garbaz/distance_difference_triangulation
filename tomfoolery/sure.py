def sqrt(x):
    pass


def triang_x(x2, x3, y3, a, b):
    q = (
        4 * x2 * x2 * a * a
        + 4 * x3 * x3 * a * a
        + 4 * y3 * y3 * a * a
        - 8 * x2 * x3 * a * a
        + 8 * b * x2 * x2 * a
        - 8 * b * x2 * x3 * a
        + 4 * b * b * x2 * x2
        - 4 * x2 * x2 * y3 * y3
    )
    return -(
        (4 * x2 * x2 * a * a * a * a) / q
        + (4 * x3 * x3 * a * a * a * a) / q
        + (4 * y3 * y3 * a * a * a * a) / q
        - (8 * x2 * x3 * a * a * a * a) / q
        + (12 * b * x2 * x2 * a * a * a) / q
        + (8 * b * x3 * x3 * a * a * a) / q
        + (8 * b * y3 * y3 * a * a * a) / q
        - (20 * b * x2 * x3 * a * a * a) / q
        - (4 * x2 * x3 * x3 * x3 * a * a) / q
        + (12 * b * b * x2 * x2 * a * a) / q
        + (8 * x2 * x2 * x3 * x3 * a * a) / q
        - (4 * x2 * x3 * y3 * y3 * a * a) / q
        - (4 * x2 * x2 * x2 * x3 * a * a) / q
        - (12 * b * b * x2 * x3 * a * a) / q
        - a * a
        - 2 * b * a
        + (4 * b * b * b * x2 * x2 * a) / q
        + (4 * b * x2 * x2 * x3 * x3 * a) / q
        - (4 * b * x2 * x2 * y3 * y3 * a) / q
        - (4 * b * x2 * x2 * x2 * x3 * a) / q
        - x2 * x2
        + (
            4
            * sqrt(
                y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                - a * a * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                - b * b * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                + x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                - 2 * a * b * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                - 2 * x3 * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                - 2 * x3 * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                + 2 * a * a * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                + 2 * b * b * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                + 4 * a * b * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                + y3 * y3 * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                - 2 * a * a * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                - 2 * b * b * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                + 2 * x3 * x3 * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                - 2 * a * b * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                + a * a * a * a * y3 * y3 * x2 * x2 * x2 * x2
                + b * b * b * b * y3 * y3 * x2 * x2 * x2 * x2
                + x3 * x3 * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                + 2 * a * b * b * b * y3 * y3 * x2 * x2 * x2 * x2
                + 2 * a * a * b * b * y3 * y3 * x2 * x2 * x2 * x2
                - 2 * a * a * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                - 2 * b * b * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                - 2 * a * b * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                + 2 * a * a * a * b * y3 * y3 * x2 * x2 * x2 * x2
                + 2 * a * a * x3 * y3 * y3 * y3 * y3 * x2 * x2 * x2
                + 2 * a * a * x3 * x3 * x3 * y3 * y3 * x2 * x2 * x2
                - 2 * a * a * a * a * x3 * y3 * y3 * x2 * x2 * x2
                - 2 * a * a * b * b * x3 * y3 * y3 * x2 * x2 * x2
                - 4 * a * a * a * b * x3 * y3 * y3 * x2 * x2 * x2
                - a * a * y3 * y3 * y3 * y3 * y3 * y3 * x2 * x2
                + a * a * a * a * y3 * y3 * y3 * y3 * x2 * x2
                + 2 * a * a * b * b * y3 * y3 * y3 * y3 * x2 * x2
                - 2 * a * a * x3 * x3 * y3 * y3 * y3 * y3 * x2 * x2
                + 2 * a * a * a * b * y3 * y3 * y3 * y3 * x2 * x2
                - a * a * b * b * b * b * y3 * y3 * x2 * x2
                - a * a * x3 * x3 * x3 * x3 * y3 * y3 * x2 * x2
                - 2 * a * a * a * b * b * b * y3 * y3 * x2 * x2
                - a * a * a * a * b * b * y3 * y3 * x2 * x2
                + a * a * a * a * x3 * x3 * y3 * y3 * x2 * x2
                + 2 * a * a * b * b * x3 * x3 * y3 * y3 * x2 * x2
                + 2 * a * a * a * b * x3 * x3 * y3 * y3 * x2 * x2
            )
            * a
        )
        / q
    ) / (2 * x2)


def triang_y(x2, x3, y3, a, b):
    q = (
        4 * x2 * x2 * a * a
        + 4 * x3 * x3 * a * a
        + 4 * y3 * y3 * a * a
        - 8 * x2 * x3 * a * a
        + 8 * b * x2 * x2 * a
        - 8 * b * x2 * x3 * a
        + 4 * b * b * x2 * x2
        - 4 * x2 * x2 * y3 * y3
    )

    return -(
        -x2 * a * a
        + x3 * a * a
        - 2 * b * x2 * a
        + 2 * b * x3 * a
        - (
            x2
            * (
                -4 * x2 * x2 * a * a * a
                - 4 * x3 * x3 * a * a * a
                - 4 * y3 * y3 * a * a * a
                + 8 * x2 * x3 * a * a * a
                - 12 * b * x2 * x2 * a * a
                - 8 * b * x3 * x3 * a * a
                - 8 * b * y3 * y3 * a * a
                + 20 * b * x2 * x3 * a * a
                + 4 * x2 * x3 * x3 * x3 * a
                - 12 * b * b * x2 * x2 * a
                - 8 * x2 * x2 * x3 * x3 * a
                + 4 * x2 * x3 * y3 * y3 * a
                + 4 * x2 * x2 * x2 * x3 * a
                + 12 * b * b * x2 * x3 * a
                - 4 * b * b * b * x2 * x2
                - 4 * b * x2 * x2 * x3 * x3
                + 4 * b * x2 * x2 * y3 * y3
                + 4 * b * x2 * x2 * x2 * x3
                - 4
                * sqrt(
                    y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                    - a * a * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                    - b * b * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                    + x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                    - 2 * a * b * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                    - 2 * x3 * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                    - 2 * x3 * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                    + 2 * a * a * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                    + 2 * b * b * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                    + 4 * a * b * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                    + y3 * y3 * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * a * a * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * b * b * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                    + 2 * x3 * x3 * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * a * b * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                    + a * a * a * a * y3 * y3 * x2 * x2 * x2 * x2
                    + b * b * b * b * y3 * y3 * x2 * x2 * x2 * x2
                    + x3 * x3 * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                    + 2 * a * b * b * b * y3 * y3 * x2 * x2 * x2 * x2
                    + 2 * a * a * b * b * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * a * a * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * b * b * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * a * b * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                    + 2 * a * a * a * b * y3 * y3 * x2 * x2 * x2 * x2
                    + 2 * a * a * x3 * y3 * y3 * y3 * y3 * x2 * x2 * x2
                    + 2 * a * a * x3 * x3 * x3 * y3 * y3 * x2 * x2 * x2
                    - 2 * a * a * a * a * x3 * y3 * y3 * x2 * x2 * x2
                    - 2 * a * a * b * b * x3 * y3 * y3 * x2 * x2 * x2
                    - 4 * a * a * a * b * x3 * y3 * y3 * x2 * x2 * x2
                    - a * a * y3 * y3 * y3 * y3 * y3 * y3 * x2 * x2
                    + a * a * a * a * y3 * y3 * y3 * y3 * x2 * x2
                    + 2 * a * a * b * b * y3 * y3 * y3 * y3 * x2 * x2
                    - 2 * a * a * x3 * x3 * y3 * y3 * y3 * y3 * x2 * x2
                    + 2 * a * a * a * b * y3 * y3 * y3 * y3 * x2 * x2
                    - a * a * b * b * b * b * y3 * y3 * x2 * x2
                    - a * a * x3 * x3 * x3 * x3 * y3 * y3 * x2 * x2
                    - 2 * a * a * a * b * b * b * y3 * y3 * x2 * x2
                    - a * a * a * a * b * b * y3 * y3 * x2 * x2
                    + a * a * a * a * x3 * x3 * y3 * y3 * x2 * x2
                    + 2 * a * a * b * b * x3 * x3 * y3 * y3 * x2 * x2
                    + 2 * a * a * a * b * x3 * x3 * y3 * y3 * x2 * x2
                )
            )
            * a
        )
        / q
        + (
            x3
            * (
                -4 * x2 * x2 * a * a * a
                - 4 * x3 * x3 * a * a * a
                - 4 * y3 * y3 * a * a * a
                + 8 * x2 * x3 * a * a * a
                - 12 * b * x2 * x2 * a * a
                - 8 * b * x3 * x3 * a * a
                - 8 * b * y3 * y3 * a * a
                + 20 * b * x2 * x3 * a * a
                + 4 * x2 * x3 * x3 * x3 * a
                - 12 * b * b * x2 * x2 * a
                - 8 * x2 * x2 * x3 * x3 * a
                + 4 * x2 * x3 * y3 * y3 * a
                + 4 * x2 * x2 * x2 * x3 * a
                + 12 * b * b * x2 * x3 * a
                - 4 * b * b * b * x2 * x2
                - 4 * b * x2 * x2 * x3 * x3
                + 4 * b * x2 * x2 * y3 * y3
                + 4 * b * x2 * x2 * x2 * x3
                - 4
                * sqrt(
                    y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                    - a * a * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                    - b * b * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                    + x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                    - 2 * a * b * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                    - 2 * x3 * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                    - 2 * x3 * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                    + 2 * a * a * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                    + 2 * b * b * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                    + 4 * a * b * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                    + y3 * y3 * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * a * a * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * b * b * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                    + 2 * x3 * x3 * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * a * b * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                    + a * a * a * a * y3 * y3 * x2 * x2 * x2 * x2
                    + b * b * b * b * y3 * y3 * x2 * x2 * x2 * x2
                    + x3 * x3 * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                    + 2 * a * b * b * b * y3 * y3 * x2 * x2 * x2 * x2
                    + 2 * a * a * b * b * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * a * a * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * b * b * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * a * b * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                    + 2 * a * a * a * b * y3 * y3 * x2 * x2 * x2 * x2
                    + 2 * a * a * x3 * y3 * y3 * y3 * y3 * x2 * x2 * x2
                    + 2 * a * a * x3 * x3 * x3 * y3 * y3 * x2 * x2 * x2
                    - 2 * a * a * a * a * x3 * y3 * y3 * x2 * x2 * x2
                    - 2 * a * a * b * b * x3 * y3 * y3 * x2 * x2 * x2
                    - 4 * a * a * a * b * x3 * y3 * y3 * x2 * x2 * x2
                    - a * a * y3 * y3 * y3 * y3 * y3 * y3 * x2 * x2
                    + a * a * a * a * y3 * y3 * y3 * y3 * x2 * x2
                    + 2 * a * a * b * b * y3 * y3 * y3 * y3 * x2 * x2
                    - 2 * a * a * x3 * x3 * y3 * y3 * y3 * y3 * x2 * x2
                    + 2 * a * a * a * b * y3 * y3 * y3 * y3 * x2 * x2
                    - a * a * b * b * b * b * y3 * y3 * x2 * x2
                    - a * a * x3 * x3 * x3 * x3 * y3 * y3 * x2 * x2
                    - 2 * a * a * a * b * b * b * y3 * y3 * x2 * x2
                    - a * a * a * a * b * b * y3 * y3 * x2 * x2
                    + a * a * a * a * x3 * x3 * y3 * y3 * x2 * x2
                    + 2 * a * a * b * b * x3 * x3 * y3 * y3 * x2 * x2
                    + 2 * a * a * a * b * x3 * x3 * y3 * y3 * x2 * x2
                )
            )
            * a
        )
        / q
        - x2 * x3 * x3
        - x2 * y3 * y3
        - b * b * x2
        + x2 * x2 * x3
        - (
            b
            * x2
            * (
                -4 * x2 * x2 * a * a * a
                - 4 * x3 * x3 * a * a * a
                - 4 * y3 * y3 * a * a * a
                + 8 * x2 * x3 * a * a * a
                - 12 * b * x2 * x2 * a * a
                - 8 * b * x3 * x3 * a * a
                - 8 * b * y3 * y3 * a * a
                + 20 * b * x2 * x3 * a * a
                + 4 * x2 * x3 * x3 * x3 * a
                - 12 * b * b * x2 * x2 * a
                - 8 * x2 * x2 * x3 * x3 * a
                + 4 * x2 * x3 * y3 * y3 * a
                + 4 * x2 * x2 * x2 * x3 * a
                + 12 * b * b * x2 * x3 * a
                - 4 * b * b * b * x2 * x2
                - 4 * b * x2 * x2 * x3 * x3
                + 4 * b * x2 * x2 * y3 * y3
                + 4 * b * x2 * x2 * x2 * x3
                - 4
                * sqrt(
                    y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                    - a * a * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                    - b * b * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                    + x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                    - 2 * a * b * y3 * y3 * x2 * x2 * x2 * x2 * x2 * x2
                    - 2 * x3 * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                    - 2 * x3 * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                    + 2 * a * a * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                    + 2 * b * b * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                    + 4 * a * b * x3 * y3 * y3 * x2 * x2 * x2 * x2 * x2
                    + y3 * y3 * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * a * a * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * b * b * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                    + 2 * x3 * x3 * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * a * b * y3 * y3 * y3 * y3 * x2 * x2 * x2 * x2
                    + a * a * a * a * y3 * y3 * x2 * x2 * x2 * x2
                    + b * b * b * b * y3 * y3 * x2 * x2 * x2 * x2
                    + x3 * x3 * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                    + 2 * a * b * b * b * y3 * y3 * x2 * x2 * x2 * x2
                    + 2 * a * a * b * b * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * a * a * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * b * b * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                    - 2 * a * b * x3 * x3 * y3 * y3 * x2 * x2 * x2 * x2
                    + 2 * a * a * a * b * y3 * y3 * x2 * x2 * x2 * x2
                    + 2 * a * a * x3 * y3 * y3 * y3 * y3 * x2 * x2 * x2
                    + 2 * a * a * x3 * x3 * x3 * y3 * y3 * x2 * x2 * x2
                    - 2 * a * a * a * a * x3 * y3 * y3 * x2 * x2 * x2
                    - 2 * a * a * b * b * x3 * y3 * y3 * x2 * x2 * x2
                    - 4 * a * a * a * b * x3 * y3 * y3 * x2 * x2 * x2
                    - a * a * y3 * y3 * y3 * y3 * y3 * y3 * x2 * x2
                    + a * a * a * a * y3 * y3 * y3 * y3 * x2 * x2
                    + 2 * a * a * b * b * y3 * y3 * y3 * y3 * x2 * x2
                    - 2 * a * a * x3 * x3 * y3 * y3 * y3 * y3 * x2 * x2
                    + 2 * a * a * a * b * y3 * y3 * y3 * y3 * x2 * x2
                    - a * a * b * b * b * b * y3 * y3 * x2 * x2
                    - a * a * x3 * x3 * x3 * x3 * y3 * y3 * x2 * x2
                    - 2 * a * a * a * b * b * b * y3 * y3 * x2 * x2
                    - a * a * a * a * b * b * y3 * y3 * x2 * x2
                    + a * a * a * a * x3 * x3 * y3 * y3 * x2 * x2
                    + 2 * a * a * b * b * x3 * x3 * y3 * y3 * x2 * x2
                    + 2 * a * a * a * b * x3 * x3 * y3 * y3 * x2 * x2
                )
            )
        )
        / q
    ) / (2 * x2 * y3)
