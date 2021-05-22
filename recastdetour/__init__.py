import sys
if sys.version_info[0] == 2:
    import Py2RecastDetour as rd
else:
    from . import Py3RecastDetour as rd


class Navmesh():
    def __init__(self):
        self._navmesh = rd.Navmesh()

    def init_by_obj(self, file_path):
        pass

    def init_by_raw(self, vertices, faces):
        pass

    def build_navmesh(self):
        pass

    def get_log(self):
        pass

    def pathfind_straight(self, start, end, vertex_mode):
        pass

    def distance_to_wall(self, point):
        pass

    def raycast(self, start, end):
        pass

    def get_settings(sself):
        pass

    def set_settings(self, settings):
        pass

    def get_partition_type(self):
        pass

    def set_partition_type(self, type):
        pass

    def get_bounding_box(self):
        pass

    def save_navmesh(self, file_path):
        pass

    def load_navmesh(self, file_path):
        pass

    def get_navmesh_trianglulation(self):
        pass
