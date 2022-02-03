def parse_label(text: str) -> Union[tuple[str, str], tuple[str, None]]:
    """ find units in parentheses """
    list_text = re.split("(.*)\((.*)\)(.*)", text)
    list_text = [text for text in list_text if text]
    if len(list_text) == 1:
        return list_text[0], None
    if len(list_text) == 2:
        return list_text[0], list_text[1]
    else:
        return text, None
