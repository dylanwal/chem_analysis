
def get_subclasses(cls, depth: int | None = 0, _count: int = 0) -> set[type]:
    subs = set(cls.__subclasses__())
    if _count < depth or depth is None:
        return subs.union(*(get_subclasses(i, depth, _count) for i in subs))

    return subs


class MixinSubClassList:
    _sub_classes: set[type] = None
    _all_sub_classes: set[type] = None

    @classmethod
    def sub_classes(cls) -> set[type]:
        if cls._sub_classes is None:
            cls._sub_classes = get_subclasses(cls, depth=0)
        return cls._sub_classes

    @classmethod
    def all_sub_classes(cls) -> set[type]:
        if cls._all_sub_classes is None:
            cls._all_sub_classes = get_subclasses(cls, depth=None)
        return get_subclasses(cls, depth=None)
