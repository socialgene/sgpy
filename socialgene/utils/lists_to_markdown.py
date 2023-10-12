"""https://gist.github.com/OsKaR31415/955b166f4a286ed427f667cb21d57bfd"""


def markdown_table_from_list(array, align: str = None):
    """
    Args:
        array: The array to make into a table. Mush be a rectangular array
               (constant width and height).
        align: The alignment of the cells : 'left', 'center' or 'right'.
    """
    # make sure every elements are strings
    array = [[str(elt) for elt in line] for line in array]
    # get the width of each column
    widths = [max(len(line[i]) for line in array) for i in range(len(array[0]))]
    # make every width at least 3 colmuns, because the separator needs it
    widths = [max(w, 3) for w in widths]
    # center text according to the widths
    array = [[elt.center(w) for elt, w in zip(line, widths)] for line in array]

    # separate the header and the body
    array_head, *array_body = array

    header = "| " + " | ".join(array_head) + " |"

    # alignment of the cells
    align = str(align).lower()  # make sure `align` is a lowercase string
    if align == "none":
        # we are just setting the position of the : in the table.
        # here there are none
        border_left = "| "
        border_center = " | "
        border_right = " |"
    elif align == "center":
        border_left = "|:"
        border_center = ":|:"
        border_right = ":|"
    elif align == "left":
        border_left = "|:"
        border_center = " |:"
        border_right = " |"
    elif align == "right":
        border_left = "| "
        border_center = ":| "
        border_right = ":|"
    else:
        raise ValueError("align must be 'left', 'right' or 'center'.")
    separator = (
        border_left + border_center.join(["-" * w for w in widths]) + border_right
    )

    # body of the table
    body = [""] * len(array)  # empty string list that we fill after
    for idx, line in enumerate(array[1:]):
        # for each line, change the body at the correct index
        body[idx] = "| " + " | ".join(line) + " |"
    body = "\n".join(body)

    return header + "\n" + separator + "\n" + body


# T = [["hi", "hello", "How do you do ?"],
#      [1, 1, 1],
#      [2, 12345678, 2],
#      ['contents', 42, 73]]

# print(make_markdown_table(T, align='left'))
