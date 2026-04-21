"""Interface with `Materials Project <https://next-gen.materialsproject.org/>`_"""

import json
import os
import pathlib
import typing

from casm.tools.shared.json_io import read_optional, read_required, safe_dump

_mp_rester = None


def get_mpr():
    global _mp_rester

    if _mp_rester is None:
        from mp_api.client import MPRester

        # Read api_key from MP_API_KEY environment variable:
        api_key = os.environ.get("MP_API_KEY")

        if api_key is None:
            raise ValueError(
                "Error in casm.tools.shared.mp_utils: "
                "MP_API_KEY environment variable not set."
            )

        _mp_rester = MPRester(api_key=api_key, use_document_model=False)

    return _mp_rester


def _get_materials(
    material_ids: str,
) -> list[dict]:
    if len(material_ids) == 0:
        return []

    # print("Fetching data from Materials Project API...")
    # print("Number of material IDs to fetch:", len(material_ids))
    # print("Material IDs:", material_ids)

    mpr = get_mpr()
    _material_docs = mpr.materials.search(material_ids=material_ids)

    from monty.json import MontyEncoder

    material_docs = []
    for _doc in _material_docs:
        doc_json = json.dumps(_doc, cls=MontyEncoder)
        doc_dict = json.loads(doc_json)
        material_docs.append(doc_dict)

    # print(
    #     "Fetched",
    #     len(material_docs),
    #     " objects from Materials Project API.",
    # )

    return material_docs


class MaterialsDocCache:
    """Cache for MaterialsDoc data from Materials Project API

    Attributes
    ----------
    current_request_index: str
        The index of the currently loaded request data, as a string.
    id_to_material: dict[str, dict]
        A mapping from material ID to dict object for the currently loaded
        request data.
    requests: dict[str, list[str]]
        A mapping from request index, as a string, to list of material IDs in that
        request.

    """

    def __init__(self, cache_dir: typing.Optional[pathlib.Path] = None):
        """

        .. rubric:: Constructor

        Parameters
        ----------
        cache_dir: Optional[pathlib.Path] = None
            The directory where cache files are stored. If None, defaults to
            ~/.casm/mp_cache/
        """
        self.current_request_index: str | None = None
        """str | None: The index of the currently loaded request data."""

        self.id_to_material: dict[str, dict] = {}
        """dict[str, dict]: A mapping from material ID to Materials document
        objects for the currently loaded request data."""

        self.requests: dict[str, list[str]] = {}
        """dict[str, list[str]]: A mapping from request index, as a string, to list of 
        material IDs in that request."""

        if cache_dir is None:
            cache_dir = pathlib.Path.home() / ".casmuser" / "mp_cache"

        self.cache_dir = cache_dir
        """pathlib.Path: The directory where cache files are stored."""

        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Load or initialize requests.json:
        requests_path = self.cache_dir / "requests.json"
        self.requests = read_optional(path=requests_path, default={})

    def switch_request(self, new_request_index: int):
        key = str(new_request_index)
        self.current_request_index = key

        data_path = self.cache_dir / f"data.{key}.json"

        if not data_path.exists():
            self.id_to_material = dict()
            if key in self.requests:
                del self.requests[key]
            safe_dump(self.requests, path=self.cache_dir / "requests.json", force=True)
            return

        material_docs = read_required(data_path)
        for doc in material_docs:
            self.id_to_material[doc["material_id"]] = doc

    def get_request_index(self, material_id: str) -> str | None:
        for key, material_ids in self.requests.items():
            if material_id in material_ids:
                return key
        return None

    def get(self, material_id: str) -> dict | None:
        key = self.get_request_index(material_id)
        if key is None:
            return None
        if key != self.current_request_index:
            self.switch_request(key)
        doc = self.id_to_material.get(material_id, None)
        return doc

    def get_missing_materials(self, material_ids: list[str]) -> list[str]:

        if len(material_ids) == 0:
            return

        missing_ids = []
        for mid in material_ids:
            if self.get_request_index(mid) is None:
                missing_ids.append(mid)

        if len(missing_ids) == 0:
            return

        next_request_index = 0
        for key in self.requests.keys():
            next_request_index = max(next_request_index, int(key) + 1)
        key = str(next_request_index)

        self.requests[key] = missing_ids
        safe_dump(self.requests, path=self.cache_dir / "requests.json", force=True)

        material_docs = _get_materials(material_ids=missing_ids)
        safe_dump(material_docs, path=self.cache_dir / f"data.{key}.json", force=True)


_material_docs_cache = MaterialsDocCache()


def get_MaterialsDocCache() -> MaterialsDocCache:
    """Get the global MaterialsDocCache instance

    Returns
    -------
    MaterialsDocCache
        The global dictCache instance.
    """
    global _material_docs_cache

    return _material_docs_cache


def get_material_docs(
    material_ids: list[str],
) -> list[dict]:
    """Get MaterialsDoc data objects from Materials Project API with caching

    Parameters
    ----------
    material_ids: list[str]
        List of material IDs to retrieve.

    Returns
    -------
    list[dict]
        List of MaterialsDoc data objects corresponding to the requested material IDs.
    """
    _material_docs_cache.get_missing_materials(material_ids=material_ids)

    material_docs = []
    for id in material_ids:
        doc = _material_docs_cache.get(id)
        if doc is None:
            raise ValueError(
                f"Error in casm.tools.shared.mp_utils.get_materials: "
                f"Failed to retrieve material ID {id}."
            )
        material_docs.append(doc)
    return material_docs
