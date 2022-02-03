import logging
import difflib


class ObjList:
    """ Object List
    Attributes
    ----------
    objs: list[Any]
        List of objs
    obj: Any
        type of object being stored
    limit: int
        maximum number of reference allowed
        default is 1000
    Methods
    -------
    add(item)
        Adds item to reference list
    remove(item)
        Removes item from reference list
    """

    def __init__(self, obj, add_objs=None, _logger: logging.Logger = None, limit: int = 1000):
        """
        Parameters
        ----------
        obj:
            the object the ref
        add_objs:
            objects will be passed to self.add()
        _logger:
            error class of parent node
        limit: int
            maximum number of reference allowed
            default is 10000
        """
        self._objs = []
        self._obj = obj
        self.limit = limit
        self.count = 0

        if _logger is None:
            self._logger = logging
        else:
            self._logger = _logger

        if add_objs is not None:
            self.add(add_objs)

    def __repr__(self):
        if self.count < 4:
            return "; ".join([repr(obj) for obj in self.objs])

        return "; ".join([repr(obj) for obj in self.objs[:2]]) + "; ..."

    def __call__(self):
        return self.objs

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.objs[item]
        elif isinstance(item, str):
            index = self._get_index_from_name(item)
            return self.objs[index]
        elif isinstance(item, slice):
            return [self.objs[i] for i in range(*item.indices(len(self.objs)))]
        else:
            mes = f"{item} not found."
            self._logger.error(mes)
            raise ValueError(mes)

    def __len__(self):
        return len(self._objs)

    def __iter__(self):
        for obj in self._objs:
            yield obj

    @property
    def objs(self):
        return self._objs

    def _get_index_from_name(self, item: str) -> int:
        """ get index from name
        Given an item name, return item index in list.
        This matching is done with difflib, so slight typing errors won't result in errors.
        Parameters
        ----------
        item: str
            name of item you are trying to find
        Returns
        -------
        index: int
            index of item in self._reference list
        Raises
        ------
        Exception
            If item name not found.
        """
        values = [i.key for i in self._objs]
        text, score = difflib.get_close_matches(item, values, n=1, cutoff=0.8)
        if score > 50:
            return values.index(text)
        else:
            mes = f"'{item}' not found."
            self._logger.error(mes)
            raise ValueError(mes)

    def add(self, objs):
        """ Add
        Adds object to reference list.
        Parameters
        ----------
        objs:
            object that you want to add to the list
        Raises
        -------
        Exception
            If invalid object is provided. An object that does not lead to a valid reference.
        """
        if not isinstance(objs, list):
            objs = [objs]

        add = []
        for obj in objs:
            if not isinstance(obj, self._obj):
                raise TypeError(f"expected: {self._obj}, received: {type(obj)} ")

            if objs in self._objs:  # check if reference s not already in list, and add it.
                self._logger.warning(f"'{obj}' already in list.")
                continue

            if hasattr(obj, "id_") and obj.id_ is None:
                obj.id_ = self.count + len(add)

            add.append(obj)

        self._objs += add
        self.count += len(add)

    def remove(self, objs):
        """ Remove
        Removes object from reference list.
        Parameters
        ----------
        objs:
            object that you want to remove to the list
        """
        if not isinstance(objs, list):
            objs = [objs]

        remove = []
        for obj in objs:
            if isinstance(obj, (str, int, slice)):
                obj = self[obj]

            if objs not in self._objs:
                self._logger.error(f"'{self._obj}'is not in list, so it can't be removed.")
                continue

            remove.append(obj)

        if not remove:
            return

        # loop through 'remove list' to remove objs
        for obj in remove:
            self._objs.remove(obj)
            self.count -= 1

    def as_dict(self) -> list:
        """ Returns list of references for serialization."""
        return [obj.as_dict() for obj in self.objs]
