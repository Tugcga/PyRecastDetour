import sys
import os
from typing import List, Tuple, Optional, Dict, Any
if sys.version_info[0] == 2:
    import Py2RecastDetour as rd
else:
    if sys.version_info[1] == 6:
        from . import Py36RecastDetour as rd
    elif sys.version_info[1] == 7:
        from . import Py37RecastDetour as rd
    elif sys.version_info[1] == 8:
        from . import Py38RecastDetour as rd
    elif sys.version_info[1] == 9:
        from . import Py39RecastDetour as rd
    else:
        from . import Py310RecastDetour as rd  # type: ignore


class Navmesh():
    def __init__(self) -> None:
        self._navmesh = rd.Navmesh()

    def init_by_obj(self, file_path: str) -> None:
        '''Initialize geomery by reading *.obj file

        Input:
            file_path - path to the file with extension *.obj
        '''
        if os.path.exists(file_path):
            extension: str = os.path.splitext(file_path)[1]
            if extension == ".obj":
                self._navmesh.init_by_obj(file_path)
            else:
                print("Fail init geometry. Only *.obj files are supported")
        else:
            print("Fail init geometry. File " + file_path + " does not exist")

    def init_by_raw(self, vertices: List[float], faces: List[int]) -> None:
        '''Initialize geometry by raw data. This data contains vertex positions and vertex indexes of polygons.

        Input:
            vertices - list of floats of the length 3x(the number of vertices) in the form [x1, y1, z1, x2, y2, z2, ...], where
                xi, yi, zi - coordinates of the i-th vertex
            faces - list of integers in the form [n1, i1, i2, ..., in1, n2, i1, i2, ..., in2, ...],
                where n1, n2, ... - the number of edges in each polygon, i1, i2, ... - indexes of polygon corners
                orientation of polygons should be in clock-wise direction

        Example: the simple plane has the following data
            [1.0, 0.0, 1.0, -1.0, 0.0, 1.0, -1.0, 0.0, -1.0, 1.0, 0.0, -1.0], [4, 0, 3, 2, 1]
        '''
        if len(vertices) % 3 == 0:
            self._navmesh.init_by_raw(vertices, faces)
        else:
            print("Fail init geometry from raw data. The number of vertices coordinates should be 3*k")

    def build_navmesh(self) -> None:
        '''Generate navmesh data. Before this method the geometry should be inited.
        '''
        self._navmesh.build_navmesh()

    def get_log(self) -> str:
        '''Return the string with inetrnal log messages.

        Output:
            string with log messages
        '''
        return self._navmesh.get_log()

    def pathfind_straight(self, start: Tuple[float, float, float], end: Tuple[float, float, float], vertex_mode: int = 0) -> Optional[List[Tuple[float, float, float]]]:
        '''Return the shortest path between start and end point inside generated navmesh.

        Input:
            start - triple of floats in the form [x, y, z]
            end - triple of floats in the form [x, y, z]
            vertex_mode - define how the result path is formed
                if vertex_mode = 0 then points adden only in path corners,
                if vertex_mode = 1 then a vertex at every polygon edge crossing where area changes is added
                if vertex_mode = 2 then vertex at every polygon edge crossing is added

        Output:
            list in the from [(x1, y1, z1), ... (xn, yn, zn)] with sequences of path points
        '''
        if len(start) == 3 and len(end) == 3:
            coordinates: List[float] = self._navmesh.pathfind_straight(start, end, vertex_mode)
            points_count = len(coordinates) // 3
            return [(coordinates[3*i], coordinates[3*i + 1], coordinates[3*i + 2]) for i in range(points_count)]
        else:
            print("Fail to find straight path. Points should be triples")
            return None

    def distance_to_wall(self, point: Tuple[float, float, float]) -> Optional[float]:
        '''Return the minimal distance between input point and navmesh edge

        Input:
            point - triple of floats in the form [x, y, z]

        Output:
            minimal distance as float number
        '''
        if len(point) == 3:
            return self._navmesh.distance_to_wall(point)
        else:
            print("Fail calculate distance to wall. The point should be triple")
            return None

    def raycast(self, start: Tuple[float, float, float], end: Tuple[float, float, float]) -> Optional[List[Tuple[float, float, float]]]:
        '''Return the segment of the line between start point and navmesh edge (or end point, if there are no collisions with navmesh edges)

        Input:
            start - triple of floats in the form [x, y, z]
            end - triple of floats in the form [x, y, z]

        Output:
            the pair [(x1, y1, z1), (x2, y2, z2)], where (x1, y1, z1) - coordinates of the start point, (x2, y2, z2) - coordinates of the finish point
        '''
        if len(start) == 3 and len(end) == 3:
            c = self._navmesh.raycast(start, end)
            if len(c) > 0:
                return [(c[0], c[1], c[2]), (c[3], c[4], c[5])]
            else:
                return None  # if calculations are fail
        else:
            print("Fails to raycast. Start and end should be triples")
            return None

    def hit_mesh(self, start: Tuple[float, float, float], end: Tuple[float, float, float]) -> Optional[Tuple[float, float, float]]:
        '''Return coordinates of the intersection point of the ray from start to end and geometry polygons

        Input:
            start - triple of floats in the form [x, y, z]
            end - triple of floats in the form [x, y, z]

        Output:
            the tuple (x, y, z) with coordinates of the intersection point
                if there are no intersections, then return coordinates of the end point
        '''
        if len(start) == 3 and len(end) == 3:
            c = self._navmesh.hit_mesh(start, end)
            if len(c) == 3:
                return (c[0], c[1], c[2])
            else:
                return None  # if calculations are fail
        else:
            print("Fails to hit mesh. Start and end should be triples")
            return None

    def get_settings(self) -> Dict[str, Any]:
        '''Return current setting, which will be used for building navmesh

        Output:
            dictionary with the following keys:
                cellSize - cell size in world units
                cellHeight - cell height in world units
                agentHeight - agent height in world units
                agentRadius - agent radius in world units
                agentMaxClimb - agent max climb in world units
                agentMaxSlope - agent max slope in degrees
                regionMinSize - region minimum size in voxels
                regionMergeSize - region merge size in voxels
                edgeMaxLen - eEdge max length in world units
                edgeMaxError - edge max error in voxels
                vertsPerPoly - the maximum number of vertices in each polygon
                detailSampleDist - detail sample distance in voxels
                detailSampleMaxError - detail sample max error in voxel heights
        '''
        return self._navmesh.get_settings()

    def set_settings(self, settings: Dict[str, Any]) -> None:
        '''Set settings for building navmesh

        Input:
            settings - dictionary with the following keys:
                cellSize - cell size in world units, should be >= 0.0001
                cellHeight - cell height in world units, should be >= 0.0001
                agentHeight - agent height in world units, should be >= 0.0
                agentRadius - agent radius in world units, should be >= 0.0
                agentMaxClimb - agent max climb in world units
                agentMaxSlope - agent max slope in degrees
                regionMinSize - region minimum size in voxels
                regionMergeSize - region merge size in voxels
                edgeMaxLen - eEdge max length in world units
                edgeMaxError - edge max error in voxels
                vertsPerPoly - the maximum number of vertices in each polygon, should be integer from 3 to 6
                detailSampleDist - detail sample distance in voxels
                detailSampleMaxError - detail sample max error in voxel heights
        '''
        self._navmesh.set_settings(settings)

    def get_partition_type(self) -> int:
        '''Retrun the index of the partition type, which used for generating polygons in the navmesh

        Output:
            one integer from 0 to 2
                0 - SAMPLE_PARTITION_WATERSHED
                1 - SAMPLE_PARTITION_MONOTONE
                2 - SAMPLE_PARTITION_LAYERS
        '''
        return self._navmesh.get_partition_type()

    def set_partition_type(self, type: int) -> None:
        '''Set partition type for generation navmesh

        Input:
            type - integer from 0 to 2
                0 - SAMPLE_PARTITION_WATERSHED
                1 - SAMPLE_PARTITION_MONOTONE
                2 - SAMPLE_PARTITION_LAYERS
        '''
        self._navmesh.set_partition_type(type)

    def get_bounding_box(self) -> Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]:
        '''Return bounding box of the mesh

        Output:
            tuple in the form (b_min, b_max), where
                b_min is a triple (x, y, z) with the lowerest corner of the bounding box
                b_max is a triple (x, y, z) with the highest corner of the bounding box
        '''
        b = self._navmesh.get_bounding_box()
        if len(b) == 6:
            return ((b[0], b[1], b[2]), (b[3], b[4], b[5]))
        else:
            return None

    def save_navmesh(self, file_path: str) -> None:
        '''Save generated navmesh to the bindary firle with extension *.bin

        Input:
            file_path - full path to the file to save
        '''
        self._navmesh.save_navmesh(file_path)

    def _generate_plane(self, plane_size: float) -> Tuple[List[float], List[int]]:
        '''Internal function, generate geometry data for simple plane
        '''
        plane_verts: List[float] = [plane_size, 0.0, plane_size, 
                                    -plane_size, 0.0, plane_size, 
                                    -plane_size, 0.0, -plane_size, 
                                    plane_size, 0.0, -plane_size]
        plane_polys: List[int] = [4, 0, 3, 2, 1]
        return (plane_verts, plane_polys)


    def load_navmesh(self, file_path: str) -> None:
        '''Load navmesh from *.bin file

        Input:
            file_path - path to the file with extension *.bin
        '''
        if os.path.exists(file_path):
            # clear generated navmesh and load simple plane
            # by default we will use the size 4.0
            self.init_by_raw(*self._generate_plane(4.0))
            settings: Dict[str, Any] = self.get_settings()
            cell_size: float = settings["cellSize"]
            if 4.0 / cell_size < 20.0:
                # regenerate the plane
                plane_size: float = cell_size * 20.0  # assume that the plane contains 20 tiles in each direction
                self.init_by_raw(*self._generate_plane(plane_size))
            # build this simple navmesh
            self.build_navmesh()
            self._navmesh.load_navmesh(file_path)
        else:
            print("Fails to load navmesh. The file " + file_path + " does not exist")

    def get_navmesh_trianglulation(self) -> Tuple[List[float], List[int]]:
        '''Return triangulation data of the generated navmesh

        Output:
            the tuple (vertices, triangles), where
                vertices is a list [x1, y1, z1, x2, y2, z2, ...] with coordinates of the navmesh vertices
                triangles is a list [t11, t12, t13, t21, t22, t23, t31, t32, t33, t41, t42, t43, ...], where ti1, ti2 and ti3 are vertex indexes of the i-th triangle
        '''
        return self._navmesh.get_navmesh_trianglulation()

    def get_navmesh_poligonization(self) -> Tuple[List[float], List[int], List[int]]:
        '''Return polygon description of the navigation mesh

        Output:
            the tuple (vertices, polygons, sizes), where
                vertices is a list [x1, y1, z1, x2, y2, z2, ...] with coordinates of the navmesh vertices
                polygons is a list [p11, p12, p13, ..., p1n1, p21, p22, p23, ..., p2n2, ...], where pij is the j-th vertex index of the polygon pi
                sizes is a list [n1, n2, ...], where ni is a size of the i-th polygon
        '''
        return self._navmesh.get_navmesh_polygonization()
