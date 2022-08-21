## What is it

This is Python binary library for generating and using navigation mesh. Based on [RecastNavigation](https://github.com/recastnavigation/recastnavigation) c++ library.

Feature:
* Generate navmesh from *.obj file and also from raw vertices and polygons data
* Save and load navmesh to (and from) binary file
* Generate triangulation and poligonization of the navmesh for using in external applications
* Find the shortest path between two input points and the shortest distance to the navmesh wall from the input point
* Works in Python 3.6 - 3.10 and Python 2.7

## How to use

Import module

```python
import recastdetour as rd
```

Create navmesh object

```python
navmesh = rd.Navmesh()
```

Init geometry

```python
navmesh.init_by_obj("location.obj")
```

Generate navmesh data

```python
navmesh.build_navmesh()
```

Find the shortest path between two points

```python
start = [2.0, 0.0, 3.0]
end = [-1.0, 0.0, -2.0]
path = navmesh.pathfind_straight(start, end)
```

For more details, see ```Navmesh``` class methods descriptions.

## Sources

Source files of the Python module is placed on the separate [repository](https://github.com/Tugcga/PyRecastDetour-Sources).