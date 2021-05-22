import sys
import os
if sys.version_info[0] == 2:
    import Py2RecastDetour as rd
else:
    from . import Py3RecastDetour as rd


class Navmesh():
    def __init__(self):
        self._navmesh = rd.Navmesh()

    def init_by_obj(self, file_path):
        if os.path.exists(file_path):
            extension = os.path.splitext(file_path)[1]
            if extension == ".obj":
                self._navmesh.init_by_obj(file_path)
            else:
                print("Fail init geometry. Only *.obj files are supported")
        else:
            print("Fail init geometry. File " + file_path + " does not exist")

    def init_by_raw(self, vertices, faces):
        if len(vertices) % 3 == 0:
            self._navmesh.init_by_raw(vertices, faces)
        else:
            print("Fail init geometry from raw data. The number of vertices coordinates should be 3*k")

    def build_navmesh(self):
        self._navmesh.build_navmesh()

    def get_log(self):
        return self._navmesh.get_log()

    def pathfind_straight(self, start, end, vertex_mode=0):
        if len(start) == 3 and len(end) == 3:
            coordinates = self._navmesh.pathfind_straight(start, end, vertex_mode)
            points_count = len(coordinates) // 3
            return [(coordinates[3*i], coordinates[3*i + 1], coordinates[3*i + 2]) for i in range(points_count)]
        else:
            print("Fail to find straight path. Points should be triples")
            return None

    def distance_to_wall(self, point):
        if len(point) == 3:
            return self._navmesh.distance_to_wall(point)
        else:
            print("Fail calculate distance to wall. The point should be triple")
            return None

    def raycast(self, start, end):
        if len(start) == 3 and len(end) == 3:
            c = self._navmesh.raycast(start, end)
            if len(c) > 0:
                return [(c[0], c[1], c[2]), (c[3], c[4], c[5])]
            else:
                return None  # if calculations are fail
        else:
            print("Fails to raycast. Start and end should be triples")
            return None

    def hit_mesh(self, start, end):
        if len(start) == 3 and len(end) == 3:
            c = self._navmesh.hit_mesh(start, end)
            if len(c) == 3:
                return tuple(c)
            else:
                return None  # if calculations are fail
        else:
            print("Fails to hit mesh. Start and end should be triples")
            return None

    def get_settings(self):
        return self._navmesh.get_settings()

    def set_settings(self, settings):
        self._navmesh.set_settings(settings)

    def get_partition_type(self):
        return self._navmesh.get_partition_type()

    def set_partition_type(self, type):
        self._navmesh.set_partition_type(type)

    def get_bounding_box(self):
        b = self._navmesh.get_bounding_box()
        if len(b) == 6:
            return ((b[0], b[1], b[2]), (b[3], b[4], b[5]))
        else:
            return None

    def save_navmesh(self, file_path):
        self._navmesh.save_navmesh(file_path)

    def load_navmesh(self, file_path):
        if os.path.exists(file_path):
            self._navmesh.load_navmesh(file_path)
        else:
            print("Fails to load navmesh. The file " + file_path + " does not exist")

    def get_navmesh_trianglulation(self):
        return self._navmesh.get_navmesh_trianglulation()
