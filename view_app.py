from PySide6 import QtWidgets, QtGui, QtCore, Qt3DRender
from PySide6.QtGui import QVector3D, QVector4D
from PySide6.Qt3DCore import Qt3DCore
from PySide6.Qt3DExtras import Qt3DExtras

import math
import struct
import os
import time
from typing import List, Tuple, Optional, Dict, Any

import recastdetour as rd

SETTINGS_SIZE: int = 256
BACKGROUND_COLOR: Tuple[int, int, int] = (55, 55, 55)
CAMERA_SPEED: float = 0.002
CAMERA_ZOOM_SPEED: float = 0.01
CAMERA_WHEEL_SPEED: float = 0.005
CAMERA_INVERT_ZOOM: bool = False
CAMERA_PAN_SPEED: float = 0.002
CAMERA_FOV: float = 45.0
CAMERA_DISTANCE_MIN: float = 0.1
CAMERA_V_MIN: float = -1.5
CAMERA_V_MAX: float = 1.5

class Lines(Qt3DRender.Qt3DRender.QGeometryRenderer):
    # each segment is an array of 3-tuples with vertex coordinates
    # so, the minimal length of each segment is 2
    def __init__(self, segments: List[List[Tuple[float, float, float]]]):
        super(Lines, self).__init__()

        points_count: int = 0  # the nmber of different points for draw
        segments_count: int = 0  # how many intervals in the lines
        for segment in segments:
            points_count += len(segment)
            segments_count += len(segment) - 1  # three points define two segments

        self._geom_view = Qt3DCore.QGeometryView(self)
        self._geom = Qt3DCore.QGeometry(self._geom_view)

        self._pos_atttr = Qt3DCore.QAttribute(self._geom)
        self._index_atttr = Qt3DCore.QAttribute(self._geom)

        self._vertex_buffer = Qt3DCore.QBuffer(self._geom)
        self._index_buffer = Qt3DCore.QBuffer(self._geom)

        self._pos_atttr.setName(Qt3DCore.QAttribute.defaultPositionAttributeName())
        self._pos_atttr.setVertexBaseType(Qt3DCore.QAttribute.Float)
        self._pos_atttr.setVertexSize(3)
        self._pos_atttr.setAttributeType(Qt3DCore.QAttribute.VertexAttribute)
        self._pos_atttr.setBuffer(self._vertex_buffer)
        self._pos_atttr.setByteStride(12)  # 3 float values (each value 4 bytes)
        self._pos_atttr.setCount(points_count)

        self._index_atttr.setAttributeType(Qt3DCore.QAttribute.IndexAttribute)
        self._index_atttr.setVertexBaseType(Qt3DCore.QAttribute.UnsignedShort)
        self._index_atttr.setBuffer(self._index_buffer)
        self._index_atttr.setCount(segments_count * 2)  # for each segment we should write two value (start and end vertex index)

        self._vertices_array = QtCore.QByteArray()
        self._indices_array = QtCore.QByteArray()
        vertices_bytes_list: List[int] = []
        indices_bytes_list: List[int] = []
        start_point_index = 0
        for segment in segments:
            # write point coordinates
            # each segment is array of 3-tuples (points)
            for point_index, point in enumerate(segment):
                for i in range(3):
                    vertices_bytes_list.extend(struct.pack("<f", point[i]))
                if point_index > 0:
                    # write previous and current point index
                    indices_bytes_list.extend(struct.pack("<H", start_point_index + point_index - 1))
                    indices_bytes_list.extend(struct.pack("<H", start_point_index + point_index))
            start_point_index += len(segment)
        self._vertices_array = QtCore.QByteArray(bytes(vertices_bytes_list))
        self._indices_array = QtCore.QByteArray(bytes(indices_bytes_list))

        # assign arrays to buffers
        self._vertex_buffer.setData(self._vertices_array)
        self._index_buffer.setData(self._indices_array)

        # add attributes to the geometry
        self._geom.addAttribute(self._pos_atttr)
        self._geom.addAttribute(self._index_atttr)

        # set geometry
        self._geom_view.setGeometry(self._geom)
        self.setView(self._geom_view)
        self._geom_view.setPrimitiveType(Qt3DCore.QGeometryView.Lines)


class PolygonMesh(Qt3DRender.Qt3DRender.QGeometryRenderer):
    def __init__(self, vertices: List[float], polygons: List[int], sizes: List[int]) -> None:
        def get_normal(i0: int, i1: int, i2: int) -> Tuple[float, float, float]:
            # return x, y, z fot the normal in the vertex i1 near vertices i0 and i2
            a: Tuple[float, float, float] = (vertices[3*i1] - vertices[3*i0], vertices[3*i1 + 1] - vertices[3*i0 + 1], vertices[3*i1 + 2] - vertices[3*i0 + 2])
            b: Tuple[float, float, float] = (vertices[3*i2] - vertices[3*i0], vertices[3*i2 + 1] - vertices[3*i0 + 1], vertices[3*i2 + 2] - vertices[3*i0 + 2])
            # get cross product
            ab: Tuple[float, float, float] = (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
            ab_length: float = math.sqrt(ab[0]**2 + ab[1]**2 + ab[2]**2)

            return (ab[0] / ab_length, ab[1] / ab_length, ab[2] / ab_length)

        super(PolygonMesh, self).__init__()
        self._geom_view = Qt3DCore.QGeometryView(self)
        self._geom = Qt3DCore.QGeometry(self._geom_view)

        self._pos_atttr = Qt3DCore.QAttribute(self._geom)
        self._normal_atttr = Qt3DCore.QAttribute(self._geom)
        self._index_atttr = Qt3DCore.QAttribute(self._geom)

        self._vertex_buffer = Qt3DCore.QBuffer(self._geom)
        self._index_buffer = Qt3DCore.QBuffer(self._geom)

        stride: int = (3 + 3) * 4  # position and normal, values are float (4 bytes)
        vertices_count: int = len(vertices) // 3
        # we should calculate the normal for each vertex
        raw_normals: List[List[float]] = [[] for i in range(vertices_count)]  # we accumulate normals for each vertex and then get average value
        triangles_count: int = 0
        triangles: List[int] = []  # store here vertex indices of the triangles
        p_index: int = 0
        for s in sizes:
            triangles_count += s - 2
            # polygon corners from p_index to p_index + s
            for i in range(1, s - 1):
                triangles.append(polygons[p_index])
                triangles.append(polygons[p_index + i])
                triangles.append(polygons[p_index + i + 1])
                raw_normals[polygons[p_index + i]].extend(get_normal(polygons[p_index + i - 1], polygons[p_index + i], polygons[p_index + i + 1]))
            # also for the first and the last vertex of the polygon
            raw_normals[polygons[p_index]].extend(get_normal(polygons[p_index + s - 1], polygons[p_index], polygons[p_index + 1]))
            raw_normals[polygons[p_index + s - 1]].extend(get_normal(polygons[p_index + s - 2], polygons[p_index + s - 1], polygons[p_index]))

            p_index += s
        normals: List[float] = []
        for i in range(vertices_count):
            v_normals: List[float] = raw_normals[i]
            vn: List[float] = [0.0, 0.0, 0.0]
            for p in range(len(v_normals) // 3):
                vn[0] += v_normals[3*p]
                vn[1] += v_normals[3*p + 1]
                vn[2] += v_normals[3*p + 2]
            vn_length: float = math.sqrt(vn[0]**2 + vn[1]**2 + vn[2]**2)
            for j in range(3):
                normals.append(vn[j] / vn_length)

        # setup attributes
        self._pos_atttr.setName(Qt3DCore.QAttribute.defaultPositionAttributeName())
        self._pos_atttr.setVertexBaseType(Qt3DCore.QAttribute.Float)
        self._pos_atttr.setVertexSize(3)
        self._pos_atttr.setAttributeType(Qt3DCore.QAttribute.VertexAttribute)
        self._pos_atttr.setBuffer(self._vertex_buffer)
        self._pos_atttr.setByteStride(stride)
        self._pos_atttr.setCount(vertices_count)

        self._normal_atttr.setName(Qt3DCore.QAttribute.defaultNormalAttributeName())
        self._normal_atttr.setVertexBaseType(Qt3DCore.QAttribute.Float)
        self._normal_atttr.setVertexSize(3)
        self._normal_atttr.setAttributeType(Qt3DCore.QAttribute.VertexAttribute)
        self._normal_atttr.setBuffer(self._vertex_buffer)
        self._normal_atttr.setByteStride(stride)
        self._normal_atttr.setByteOffset(3 * 4)  # normals store after 3 values of the position (each 4 bytes)
        self._normal_atttr.setCount(vertices_count)

        self._index_atttr.setAttributeType(Qt3DCore.QAttribute.IndexAttribute)
        self._index_atttr.setVertexBaseType(Qt3DCore.QAttribute.UnsignedShort)
        self._index_atttr.setBuffer(self._index_buffer)

        self._index_atttr.setCount(triangles_count * 3)

        # next setup buffers
        # vertices
        vertices_bytes_list: List[int] = []
        for v_index in range(vertices_count):
            vertices_bytes_list.extend(struct.pack("<f", vertices[3*v_index]))
            vertices_bytes_list.extend(struct.pack("<f", vertices[3*v_index + 1]))
            vertices_bytes_list.extend(struct.pack("<f", vertices[3*v_index + 2]))

            vertices_bytes_list.extend(struct.pack("<f", normals[3*v_index]))
            vertices_bytes_list.extend(struct.pack("<f", normals[3*v_index + 1]))
            vertices_bytes_list.extend(struct.pack("<f", normals[3*v_index + 2]))
        self._vertices_array = QtCore.QByteArray(bytes(vertices_bytes_list))
        
        # indices
        indices_bytes_list: List[int] = []
        for t in triangles:
            indices_bytes_list.extend(struct.pack("<H", t))
        self._indices_array = QtCore.QByteArray(bytes(indices_bytes_list))
        
        # assign arrays to buffers
        self._vertex_buffer.setData(self._vertices_array)
        self._index_buffer.setData(self._indices_array)

        # add attributes to the geometry
        self._geom.addAttribute(self._pos_atttr)
        self._geom.addAttribute(self._normal_atttr)
        self._geom.addAttribute(self._index_atttr)

        # set geometry
        self._geom_view.setGeometry(self._geom)
        self.setView(self._geom_view)


class WireframeMaterial(Qt3DRender.Qt3DRender.QMaterial):
    def __init__(self, root_entity: Qt3DCore.QEntity) -> None:
        super(WireframeMaterial, self).__init__(root_entity)

        self._vertex_effect = Qt3DRender.Qt3DRender.QEffect()
        self._technique = Qt3DRender.Qt3DRender.QTechnique()
        self._render_pass = Qt3DRender.Qt3DRender.QRenderPass()
        self._shader = Qt3DRender.Qt3DRender.QShaderProgram()
        self._filter = Qt3DRender.Qt3DRender.QFilterKey()

        self._shader.setVertexShaderCode(QtCore.QByteArray(b'''
#version 330 core

in vec3 vertexPosition;
in vec3 vertexNormal;

out EyeSpaceVertex {
    vec3 position;
    vec3 normal;
} vs_out;

uniform mat4 modelView;
uniform mat3 modelViewNormal;
uniform mat4 mvp;

void main()
{
    vs_out.normal = normalize( modelViewNormal * vertexNormal );
    vs_out.position = vec3( modelView * vec4( vertexPosition, 1.0 ) );

    gl_Position = mvp * vec4( vertexPosition, 1.0 );
}

'''))
        self._shader.setFragmentShaderCode(QtCore.QByteArray(b'''
#version 330 core

uniform struct LightInfo {
    vec4 position;
    vec3 intensity;
} light;

uniform struct LineInfo {
    float width;
    vec4 color;
} line;

uniform vec3 ka;            // Ambient reflectivity
uniform vec3 kd;            // Diffuse reflectivity
uniform vec3 ks;            // Specular reflectivity
uniform float shininess;    // Specular shininess factor

in WireframeVertex {
    vec3 position;
    vec3 normal;
    noperspective vec4 edgeA;
    noperspective vec4 edgeB;
    flat int configuration;
} fs_in;

out vec4 fragColor;

vec3 adsModel( const in vec3 pos, const in vec3 n )
{
    // Calculate the vector from the light to the fragment
    vec3 s = normalize( vec3( light.position ) - pos );

    // Calculate the vector from the fragment to the eye position (the
    // origin since this is in "eye" or "camera" space
    vec3 v = normalize( -pos );

    // Refleft the light beam using the normal at this fragment
    vec3 r = reflect( -s, n );

    // Calculate the diffus component
    vec3 diffuse = vec3( max( dot( s, n ), 0.0 ) );

    // Calculate the specular component
    vec3 specular = vec3( pow( max( dot( r, v ), 0.0 ), shininess ) );

    // Combine the ambient, diffuse and specular contributions
    return light.intensity * ( ka + kd * diffuse + ks * specular );
}

vec4 shadeLine( const in vec4 color )
{
    // Find the smallest distance between the fragment and a triangle edge
    float d;
    if ( fs_in.configuration == 0 )
    {
        // Common configuration
        d = min( fs_in.edgeA.x, fs_in.edgeA.y );
        d = min( d, fs_in.edgeA.z );
    }
    else
    {
        // Handle configuration where screen space projection breaks down
        // Compute and compare the squared distances
        vec2 AF = gl_FragCoord.xy - fs_in.edgeA.xy;
        float sqAF = dot( AF, AF );
        float AFcosA = dot( AF, fs_in.edgeA.zw );
        d = abs( sqAF - AFcosA * AFcosA );

        vec2 BF = gl_FragCoord.xy - fs_in.edgeB.xy;
        float sqBF = dot( BF, BF );
        float BFcosB = dot( BF, fs_in.edgeB.zw );
        d = min( d, abs( sqBF - BFcosB * BFcosB ) );

        // Only need to care about the 3rd edge for some configurations.
        if ( fs_in.configuration == 1 || fs_in.configuration == 2 || fs_in.configuration == 4 )
        {
            float AFcosA0 = dot( AF, normalize( fs_in.edgeB.xy - fs_in.edgeA.xy ) );
            d = min( d, abs( sqAF - AFcosA0 * AFcosA0 ) );
        }

        d = sqrt( d );
    }

    // Blend between line color and phong color
    float mixVal;
    if ( d < line.width - 1.0 )
    {
        mixVal = 1.0;
    }
    else if ( d > line.width + 1.0 )
    {
        mixVal = 0.0;
    }
    else
    {
        float x = d - ( line.width - 1.0 );
        mixVal = exp2( -2.0 * ( x * x ) );
    }

    return mix( color, line.color, mixVal );
}

void main()
{
    // Calculate the color from the phong model
    vec4 color = vec4( adsModel( fs_in.position, normalize( fs_in.normal ) ), 1.0 );
    fragColor = shadeLine( color );
}

'''))

        self._shader.setGeometryShaderCode(QtCore.QByteArray(b'''
#version 330 core

layout( triangles ) in;
layout( triangle_strip, max_vertices = 3 ) out;

in EyeSpaceVertex {
    vec3 position;
    vec3 normal;
} gs_in[];

out WireframeVertex {
    vec3 position;
    vec3 normal;
    noperspective vec4 edgeA;
    noperspective vec4 edgeB;
    flat int configuration;
} gs_out;

uniform mat4 viewportMatrix;

const int infoA[]  = int[]( 0, 0, 0, 0, 1, 1, 2 );
const int infoB[]  = int[]( 1, 1, 2, 0, 2, 1, 2 );
const int infoAd[] = int[]( 2, 2, 1, 1, 0, 0, 0 );
const int infoBd[] = int[]( 2, 2, 1, 2, 0, 2, 1 );

vec2 transformToViewport( const in vec4 p )
{
    return vec2( viewportMatrix * ( p / p.w ) );
}

void main()
{
    gs_out.configuration = int(gl_in[0].gl_Position.z < 0) * int(4)
           + int(gl_in[1].gl_Position.z < 0) * int(2)
           + int(gl_in[2].gl_Position.z < 0);

    // If all vertices are behind us, cull the primitive
    if (gs_out.configuration == 7)
        return;

    // Transform each vertex into viewport space
    vec2 p[3];
    p[0] = transformToViewport( gl_in[0].gl_Position );
    p[1] = transformToViewport( gl_in[1].gl_Position );
    p[2] = transformToViewport( gl_in[2].gl_Position );

    if (gs_out.configuration == 0)
    {
        // Common configuration where all vertices are within the viewport
        gs_out.edgeA = vec4(0.0);
        gs_out.edgeB = vec4(0.0);

        // Calculate lengths of 3 edges of triangle
        float a = length( p[1] - p[2] );
        float b = length( p[2] - p[0] );
        float c = length( p[1] - p[0] );

        // Calculate internal angles using the cosine rule
        float alpha = acos( ( b * b + c * c - a * a ) / ( 2.0 * b * c ) );
        float beta = acos( ( a * a + c * c - b * b ) / ( 2.0 * a * c ) );

        // Calculate the perpendicular distance of each vertex from the opposing edge
        float ha = abs( c * sin( beta ) );
        float hb = abs( c * sin( alpha ) );
        float hc = abs( b * sin( alpha ) );

        // Now add this perpendicular distance as a per-vertex property in addition to
        // the position and normal calculated in the vertex shader.

        // Vertex 0 (a)
        gs_out.edgeA = vec4( ha, 0.0, 0.0, 0.0 );
        gs_out.normal = gs_in[0].normal;
        gs_out.position = gs_in[0].position;
        gl_Position = gl_in[0].gl_Position;
        EmitVertex();

        // Vertex 1 (b)
        gs_out.edgeA = vec4( 0.0, hb, 0.0, 0.0 );
        gs_out.normal = gs_in[1].normal;
        gs_out.position = gs_in[1].position;
        gl_Position = gl_in[1].gl_Position;
        EmitVertex();

        // Vertex 2 (c)
        gs_out.edgeA = vec4( 0.0, 0.0, hc, 0.0 );
        gs_out.normal = gs_in[2].normal;
        gs_out.position = gs_in[2].position;
        gl_Position = gl_in[2].gl_Position;
        EmitVertex();

        // Finish the primitive off
        EndPrimitive();
    }
    else
    {
        // Viewport projection breaks down for one or two vertices.
        // Caclulate what we can here and defer rest to fragment shader.
        // Since this is coherent for the entire primitive the conditional
        // in the fragment shader is still cheap as all concurrent
        // fragment shader invocations will take the same code path.

        // Copy across the viewport-space points for the (up to) two vertices
        // in the viewport
        gs_out.edgeA.xy = p[infoA[gs_out.configuration]];
        gs_out.edgeB.xy = p[infoB[gs_out.configuration]];

        // Copy across the viewport-space edge vectors for the (up to) two vertices
        // in the viewport
        gs_out.edgeA.zw = normalize( gs_out.edgeA.xy - p[ infoAd[gs_out.configuration] ] );
        gs_out.edgeB.zw = normalize( gs_out.edgeB.xy - p[ infoBd[gs_out.configuration] ] );

        // Pass through the other vertex attributes
        gs_out.normal = gs_in[0].normal;
        gs_out.position = gs_in[0].position;
        gl_Position = gl_in[0].gl_Position;
        EmitVertex();

        gs_out.normal = gs_in[1].normal;
        gs_out.position = gs_in[1].position;
        gl_Position = gl_in[1].gl_Position;
        EmitVertex();

        gs_out.normal = gs_in[2].normal;
        gs_out.position = gs_in[2].position;
        gl_Position = gl_in[2].gl_Position;
        EmitVertex();

        // Finish the primitive off
        EndPrimitive();
    }
}

'''))
        
        self._technique.graphicsApiFilter().setApi(Qt3DRender.Qt3DRender.QGraphicsApiFilter.OpenGL)
        self._technique.graphicsApiFilter().setMajorVersion(3)
        self._technique.graphicsApiFilter().setMinorVersion(1)
        self._technique.graphicsApiFilter().setProfile(Qt3DRender.Qt3DRender.QGraphicsApiFilter.CoreProfile)

        self._filter.setParent(self);
        self._filter.setName("renderingStyle")
        self._filter.setValue("forward")

        self._technique.addFilterKey(self._filter)
        self.line_width = Qt3DRender.Qt3DRender.QParameter("line.width", 1.0)
        self.line_color = Qt3DRender.Qt3DRender.QParameter("line.color", QVector4D(1.0, 1.0, 1.0, 1.0))
        self.light_postion = Qt3DRender.Qt3DRender.QParameter("light.position", QVector4D(0.0, 0.0, 0.0, 1.0))
        self.light_intensity = Qt3DRender.Qt3DRender.QParameter("light.intensity", QVector3D(1.0, 1.0, 1.0))
        self._technique.addParameter(self.line_width)
        self._technique.addParameter(self.line_color)
        self._technique.addParameter(self.light_postion)
        self._technique.addParameter(self.light_intensity)

        self._render_pass.setShaderProgram(self._shader)

        self._technique.addRenderPass(self._render_pass)
        self._vertex_effect.addTechnique(self._technique)

        self.setEffect(self._vertex_effect)

        self.ka_param = Qt3DRender.Qt3DRender.QParameter("ka", QVector4D(0.1, 0.1, 0.1, 1.0))
        self.kd_param = Qt3DRender.Qt3DRender.QParameter("kd", QVector4D(0.7, 0.7, 0.7, 1.0))
        self.ks_param = Qt3DRender.Qt3DRender.QParameter("ks", QVector3D(0.95, 0.95, 0.95))
        self.shininess_param = Qt3DRender.Qt3DRender.QParameter("shininess", 150.0)
        self.addParameter(self.ka_param)
        self.addParameter(self.kd_param)
        self.addParameter(self.ks_param)
        self.addParameter(self.shininess_param)

    def set_light_position(self, x: float, y: float, z: float) -> None:
        self.light_postion.setValue(QVector4D(x, y, z, 1.0))

    def set_light_intensity(self, value: float) -> None:
        self.light_intensity.setValue(QVector3D(value, value, value))

    def set_line_width(self, value: float) -> None:
        self.line_width.setValue(value)

    def set_line_color(self, r: float, g: float, b: float) -> None:
        self.line_color.setValue(QVector4D(r, g, b, 1.0))

    def set_material_ambient(self, r: float, g: float, b: float) -> None:
        self.ka_param.setValue(QVector4D(r, g, b, 1.0))

    def set_material_diffuse(self, r: float, g: float, b: float) -> None:
        self.kd_param.setValue(QVector4D(r, g, b, 1.0))

    def set_material_specular(self, r: float, g: float, b: float) -> None:
        self.ks_param.setValue(QVector3D(r, g, b))

    def set_material_shiness(self, value: float) -> None:
        self.shininess_param.setValue(value)


class ConstantAlphaMaterial(Qt3DRender.Qt3DRender.QMaterial):
    def __init__(self, root_entity: Qt3DCore.QEntity) -> None:
        super(ConstantAlphaMaterial, self).__init__(root_entity)

        self._vertex_effect = Qt3DRender.Qt3DRender.QEffect()
        self._technique = Qt3DRender.Qt3DRender.QTechnique()
        self._render_pass = Qt3DRender.Qt3DRender.QRenderPass()
        self._shader = Qt3DRender.Qt3DRender.QShaderProgram()
        self._filter = Qt3DRender.Qt3DRender.QFilterKey()

        self._shader.setVertexShaderCode(QtCore.QByteArray(b'''
#version 330 core

in vec3 vertexPosition;
uniform mat4 mvp;

void main()
{
    gl_Position = mvp * vec4( vertexPosition, 1.0 );
}

'''))
        self._shader.setFragmentShaderCode(QtCore.QByteArray(b'''
#version 330 core

uniform vec4 color;
out vec4 fragColor;

void main()
{
    fragColor = color;
}

'''))

        self._technique.graphicsApiFilter().setApi(Qt3DRender.Qt3DRender.QGraphicsApiFilter.OpenGL)
        self._technique.graphicsApiFilter().setMajorVersion(3)
        self._technique.graphicsApiFilter().setMinorVersion(1)
        self._technique.graphicsApiFilter().setProfile(Qt3DRender.Qt3DRender.QGraphicsApiFilter.CoreProfile)

        self._filter.setParent(self);
        self._filter.setName("renderingStyle")
        self._filter.setValue("forward")

        self._technique.addFilterKey(self._filter)
        self._render_pass.setShaderProgram(self._shader)

        self.depth_test = Qt3DRender.Qt3DRender.QDepthTest()
        self.blend_arguments = Qt3DRender.Qt3DRender.QBlendEquationArguments()
        self.blend_equation = Qt3DRender.Qt3DRender.QBlendEquation()

        self.depth_test.setDepthFunction(Qt3DRender.Qt3DRender.QDepthTest.Always)

        self.blend_arguments.setSourceRgb(Qt3DRender.Qt3DRender.QBlendEquationArguments.SourceAlpha)
        self.blend_arguments.setDestinationRgb(Qt3DRender.Qt3DRender.QBlendEquationArguments.OneMinusSourceAlpha)
        self.blend_equation.setBlendFunction(Qt3DRender.Qt3DRender.QBlendEquation.Add)

        self._render_pass.addRenderState(self.depth_test)
        self._render_pass.addRenderState(self.blend_arguments)
        self._render_pass.addRenderState(self.blend_equation)

        self._technique.addRenderPass(self._render_pass)
        self._vertex_effect.addTechnique(self._technique)

        self.color_param = Qt3DRender.Qt3DRender.QParameter("color", QVector4D(0.7, 0.7, 0.7, 1.0))
        self.addParameter(self.color_param)

        self.setEffect(self._vertex_effect)

    def set_color(self, r: float, g: float, b: float, a: float) -> None:
        self.color_param.setValue(QVector4D(r, g, b, a))


class CanvasWidget(Qt3DExtras.Qt3DWindow):
    def __init__(self, main) -> None:
        super(CanvasWidget, self).__init__()
        self._main = main
        frame_graph = self.activeFrameGraph()
        frame_graph.setClearColor(QtGui.QColor(*BACKGROUND_COLOR))

        self.root_entity = Qt3DCore.QEntity()

        self._geometry_material: WireframeMaterial = WireframeMaterial(self.root_entity)
        self._geometry_material.set_material_specular(0.0, 0.0, 0.0)
        self._geometry_material.set_line_color(0.8, 0.8, 0.8)
        self._geometry_material.set_line_width(0.5)
        self._geometry_material.set_material_diffuse(1.0, 1.0, 1.0)
        self._geometry_material.set_material_ambient(1.0, 1.0, 1.0)
        self._geometry_material.set_light_intensity(0.3)

        self._navmesh_material: ConstantAlphaMaterial = ConstantAlphaMaterial(self.root_entity)
        self._navmesh_material.set_color(0.3, 0.73, 0.96, 0.7)

        self.setRootEntity(self.root_entity)

        # parameters for camera controll
        # define pan and zoom sped with respect to the size of the level bouding box
        self._camera_pan_speed: float = 0.0
        self._camera_zoom_speed: float = 0.0
        self._is_left_press: bool = False
        self._is_middle_press: bool = False
        self._is_right_press: bool = False
        self._left_press_x: int = 0
        self._left_press_y: int = 0
        self._middle_press_x: int = 0
        self._middle_press_y: int = 0
        self._right_press_x: int = 0
        self._right_press_y: int = 0
        self._distance_ref: float = 16.0
        self._distance: float = self._distance_ref
        self._uv_ref: Tuple[float, float] = (math.pi / 3, math.pi / 6)
        self._uv: Tuple[float, float] = (self._uv_ref[0], self._uv_ref[1])
        self._camera_center_ref: Tuple[float, float, float] = (0.0, 0.0, 0.0)
        self._camera_center = QVector3D(*self._camera_center_ref)
        self._camera_position = QVector3D(0.0, 0.0, 0.0)
        self.camera().setPosition(self._camera_position)
        self.camera().setViewCenter(self._camera_center)
        self.camera().setUpVector(QVector3D(0, 1, 0))
        self.camera().lens().setPerspectiveProjection(CAMERA_FOV, 1.0, 0.1, 1000)
        
        self._set_camera_position()

        self._level_vertices: List[Tuple[float, float, float]] = []
        self._level: Optional[Qt3DCore.QEntity] = None
        self._navmesh: Optional[Qt3DCore.QEntity] = None
        self._navmesh_lines: Optional[Qt3DCore.QEntity] = None

    def _uv_to_xyz(self, u: float, v: float) -> Tuple[float, float, float]:
        x = math.cos(u) * math.cos(v) * self._distance
        y = math.sin(v) * self._distance
        z = math.sin(u) * math.cos(v) * self._distance

        return x, y, z

    def _set_camera_position(self) -> None:
        x, y, z = self._uv_to_xyz(*self._uv)

        self._camera_position.setX(x + self._camera_center.x())
        self._camera_position.setY(y + self._camera_center.y())
        self._camera_position.setZ(z + self._camera_center.z())
        self.camera().setViewCenter(self._camera_center)
        self.camera().setPosition(self._camera_position)

    def _get_current_uv(self, event_point: QtCore.QPointF) -> Tuple[float, float]:
        p_x, p_y = event_point.x(), event_point.y()
        x_delta: float = p_x - self._left_press_x
        y_delta: float = p_y - self._left_press_y
        delta_u: float = x_delta * math.pi * CAMERA_SPEED
        delta_v: float = y_delta * math.pi * CAMERA_SPEED
        u = self._uv_ref[0] + delta_u
        v = self._uv_ref[1] + delta_v

        if v > CAMERA_V_MAX:
            v = CAMERA_V_MAX
        elif v < CAMERA_V_MIN:
            v = CAMERA_V_MIN
        return u, v

    def _get_current_distance(self, event_point: QtCore.QPointF) -> float:
        x_delta: float = event_point.x() - self._right_press_x
        y_delta: float = event_point.y() - self._right_press_y
        d: float = self._distance_ref + y_delta * CAMERA_ZOOM_SPEED * self._camera_zoom_speed * (1.0 if CAMERA_INVERT_ZOOM else -1.0)
        if d < CAMERA_DISTANCE_MIN:
            d = CAMERA_DISTANCE_MIN
        return d

    def _get_current_center(self, event_point: QtCore.QPointF) -> Tuple[float, float, float]:
        x_delta: float = event_point.x() - self._middle_press_x
        y_delta: float = event_point.y() - self._middle_press_y
        # get camera view vector
        view_vector: Tuple[float, float, float] = self._uv_to_xyz(*self._uv)
        length: float = math.sqrt(view_vector[0]**2 + view_vector[1]**2 + view_vector[2]**2)
        view_vector = (-1 * view_vector[0] / length, -1 * view_vector[1] / length, -1 * view_vector[2] / length)
        # calculate screen x-vector as cross producе of up-vector and view-vector
        x_vector: Tuple[float, float, float] = (-1 * view_vector[2], 0.0, view_vector[0])
        # normalize it
        x_length: float = math.sqrt(x_vector[0]**2 + x_vector[1]**2 + x_vector[2]**2)
        x_vector = (x_vector[0] / x_length, x_vector[1] / x_length, x_vector[2] / x_length)
        # next y_vector is cross product of x_vector and view_vector
        y_vector: Tuple[float, float, float] = (x_vector[1]*view_vector[2] - x_vector[2]*view_vector[1], 
                    x_vector[2]*view_vector[0] - x_vector[0]*view_vector[2],
                    x_vector[0]*view_vector[1] - x_vector[1]*view_vector[0])
        shift_vector: Tuple[float, float, float] = (CAMERA_PAN_SPEED * self._camera_pan_speed * (-x_delta * x_vector[0] + y_delta * y_vector[0]),
                                                    CAMERA_PAN_SPEED * self._camera_pan_speed * (-x_delta * x_vector[1] + y_delta * y_vector[1]),
                                                    CAMERA_PAN_SPEED * self._camera_pan_speed * (-x_delta * x_vector[2] + y_delta * y_vector[2]))
        return (self._camera_center_ref[0] + shift_vector[0],
                self._camera_center_ref[1] + shift_vector[1],
                self._camera_center_ref[2] + shift_vector[2])

    def _frame_level(self) -> None:
        if len(self._level_vertices) > 0:
            min_x = float("inf")
            min_y = float("inf")
            min_z = float("inf")
            max_x = -float("inf")
            max_y = -float("inf")
            max_z = -float("inf")
            for v in self._level_vertices:
                if v[0] < min_x:
                    min_x = v[0]
                if v[0] > max_x:
                    max_x = v[0]
                if v[1] < min_y:
                    min_y = v[1]
                if v[1] > max_y:
                    max_y = v[1]
                if v[2] < min_z:
                    min_z = v[2]
                if v[2] > max_z:
                    max_z = v[2]
            self._camera_center_ref = ((min_x + max_x) / 2.0, (min_y + max_y) / 2.0, (min_z + max_z) / 2.0)
            self._camera_center.setX(self._camera_center_ref[0])
            self._camera_center.setY(self._camera_center_ref[1])
            self._camera_center.setZ(self._camera_center_ref[2])

            # also we should calculate the distance to observe all point of the bounding box
            size: float = max(max_x - min_x, max_y - min_y)
            size = max(size, max_z - min_z)
            self._distance = size / math.tan(CAMERA_FOV * math.pi / 360.0)
            self._distance_ref = self._distance
            self._set_camera_position()

            # set camera speed values
            self._camera_zoom_speed = size
            self._camera_pan_speed = size

    def _load_level(self) -> None:
        path: str = QtCore.QUrl.fromLocalFile(self._level_obj_path)
        self._level = Qt3DCore.QEntity(self.root_entity)
        self._level_mesh = Qt3DRender.Qt3DRender.QMesh()
        self._level_mesh.setSource(path)

        self._level_tfm = Qt3DCore.QTransform()
        self._level.addComponent(self._level_mesh)
        self._level.addComponent(self._geometry_material)
        self._level.addComponent(self._level_tfm)

        # parse input obj and extract vertices
        self._level_vertices.clear()
        with open(self._level_obj_path, "r") as file:
            for line in file:
                if line[0] == "v" and line[1] == " ":
                    parts: List[str] = line[2:].split(" ")
                    if len(parts) == 3:
                        self._level_vertices.append((float(parts[0]), float(parts[1]), float(parts[2])))
        self._frame_level()

    def _clear_level(self) -> None:
        if self._level is not None:
            self._level.setParent(None)  # type: ignore
        self._level_define = False
        self._level_obj_path = ""

    def _clear_navmesh(self) -> None:
        if self._navmesh is not None:
            self._navmesh.setParent(None)  # type: ignore
            self._navmesh = None
        if self._navmesh_lines is not None:
            self._navmesh_lines.setParent(None)  # type: ignore
            self._navmesh_lines = None

    def keyPressEvent(self, event: QtGui.QKeyEvent) -> None:
        if event.isAutoRepeat() is False and event.key() == 70:  # 70 = f
            self._frame_level()

    def mouseMoveEvent(self, event: QtGui.QMouseEvent) -> None:
        if self._is_left_press:
            current_point = event.position()
            self._uv = self._get_current_uv(current_point)
            self._set_camera_position()
        if self._is_middle_press:
            current_point = event.position()
            x, y, z = self._get_current_center(current_point)
            self._camera_center.setX(x)
            self._camera_center.setY(y)
            self._camera_center.setZ(z)
            self._set_camera_position()
        if self._is_right_press:
            current_point = event.position()
            self._distance = self._get_current_distance(current_point)
            self._set_camera_position()

    # events
    def mousePressEvent(self, event: QtGui.QMouseEvent) -> None:
        press_point: QtCore.QPointF
        if event.button() == QtCore.Qt.MouseButton.LeftButton:
            self._is_left_press = True
            press_point = event.position()
            self._left_press_x, self._left_press_y = int(press_point.x()), int(press_point.y())
        elif event.button() == QtCore.Qt.MouseButton.MiddleButton:
            self._is_middle_press = True
            press_point= event.position()
            self._middle_press_x, self._middle_press_y = int(press_point.x()), int(press_point.y())
        elif event.button() == QtCore.Qt.MouseButton.RightButton:
            self._is_right_press = True
            press_point = event.position()
            self._right_press_x, self._right_press_y = int(press_point.x()), int(press_point.y())

    def mouseReleaseEvent(self, event: QtGui.QMouseEvent) -> None:
        current_point: QtCore.QPointF
        if event.button() == QtCore.Qt.MouseButton.LeftButton and self._is_left_press :
            self._is_left_press = False
            current_point = event.position()
            self._uv_ref = self._get_current_uv(current_point)
            self._uv = (self._uv_ref[0], self._uv_ref[1])
            self._set_camera_position()
        if event.button() == QtCore.Qt.MouseButton.MiddleButton and self._is_middle_press:
            self._is_middle_press = False
            current_point = event.position()
            self._camera_center_ref = self._get_current_center(current_point)
        if event.button() == QtCore.Qt.MouseButton.RightButton and self._is_right_press:
            self._is_right_press = False
            current_point = event.position()
            self._distance_ref = self._get_current_distance(current_point)

    def wheelEvent(self, event: QtGui.QMouseEvent) -> None:
        delta: float = event.angleDelta().y()
        self._distance_ref += -1 * delta * CAMERA_WHEEL_SPEED * self._camera_zoom_speed * (1.0 if CAMERA_INVERT_ZOOM else -1.0)
        if self._distance_ref < CAMERA_DISTANCE_MIN:
            self._distance_ref = CAMERA_DISTANCE_MIN
        self._distance = self._distance_ref
        self._set_camera_position()

    # public methods
    def clear_scene(self) -> None:
        self._clear_navmesh()
        self._clear_level()

    def load_level_command(self, obj_path: str) -> None:
        self._clear_level()
        self._level_obj_path = obj_path
        self._level_define = True
        
        self._load_level()
        self._clear_navmesh()

    def draw_navmesh(self, vertices: List[float], polygons: List[int], sizes: List[int], is_frame: bool=False) -> None:
        self._clear_navmesh()

        if len(vertices) > 0 and len(polygons) > 0 and len(sizes)> 0:
            self._navmesh_geometry = PolygonMesh(vertices, polygons, sizes)
            self._navmesh = Qt3DCore.QEntity(self.root_entity)
            self._navmesh_tfm = Qt3DCore.QTransform()

            self._navmesh.addComponent(self._navmesh_geometry)
            self._navmesh.addComponent(self._navmesh_tfm)
            self._navmesh.addComponent(self._navmesh_material)

            self._navmesh_lines = Qt3DCore.QEntity(self.root_entity)
            self._navmesh_lines_tfm = Qt3DCore.QTransform()

            # next we should create lines for navmesh polygons
            polygon_lines: List[List[Tuple[float, float, float]]] = []
            p_index: int = 0
            for p_size in sizes:
                polygon_border: List[Tuple[float, float, float]] = []
                for i in range(p_size):
                    p: int = polygons[p_index + i]
                    polygon_border.append((vertices[3*p], vertices[3*p + 1], vertices[3*p + 2]))
                # and also add the last segment
                polygon_border.append(polygon_border[0])
                polygon_lines.append(polygon_border)
                p_index += p_size

            self._navmesh_lines_geometry = Lines(polygon_lines)
            self._navmesh_lines.addComponent(self._navmesh_lines_geometry)
            self._navmesh_lines.addComponent(self._navmesh_lines_tfm)
            self._navmesh_lines.addComponent(self._navmesh_material)

        if is_frame:
            self._level_vertices.clear()
            vertices_count: int = len(vertices) // 3
            for i in range(vertices_count):
                self._level_vertices.append((vertices[3*i], vertices[3*i + 1], vertices[3*i + 2]))
            self._frame_level()


class SettingsWidget(QtWidgets.QWidget):
    def __init__(self, canvas, main):
        super(SettingsWidget, self).__init__()
        self._navmesh: Optional[rd.Navmesh] = None
        self._navmesh_mode: int = -1  # -1 - unknown, 0 - level, 1 - binary data

        self._canvas: CanvasWidget = canvas
        self._main: ViewMainWidget = main
        self.setMinimumWidth(SETTINGS_SIZE)
        self.setMaximumWidth(SETTINGS_SIZE)

        self._regenerate_on_change: bool = False

        self._settings_cell_size = QtWidgets.QDoubleSpinBox()
        self._settings_cell_size.setValue(0.3)
        self._settings_cell_size.setSingleStep(0.05)
        self._settings_cell_size.setMinimum(0.0001)
        self._settings_cell_size.valueChanged.connect(self._auto_bake)

        self._settings_cell_height = QtWidgets.QDoubleSpinBox()
        self._settings_cell_height.setValue(0.2)
        self._settings_cell_height.setSingleStep(0.05)
        self._settings_cell_height.setMinimum(0.0001)
        self._settings_cell_height.valueChanged.connect(self._auto_bake)

        self._settings_agent_height = QtWidgets.QDoubleSpinBox()
        self._settings_agent_height.setValue(2.0)
        self._settings_agent_height.setSingleStep(0.5)
        self._settings_agent_height.setMinimum(0.0001)
        self._settings_agent_height.valueChanged.connect(self._auto_bake)

        self._settings_agent_radius = QtWidgets.QDoubleSpinBox()
        self._settings_agent_radius.setValue(0.6)
        self._settings_agent_radius.setSingleStep(0.1)
        self._settings_agent_radius.setMinimum(0.0001)
        self._settings_agent_radius.valueChanged.connect(self._auto_bake)

        self._settings_agent_max_climb = QtWidgets.QDoubleSpinBox()
        self._settings_agent_max_climb.setValue(0.9)
        self._settings_agent_max_climb.setSingleStep(0.1)
        self._settings_agent_max_climb.setMinimum(0.0001)
        self._settings_agent_max_climb.valueChanged.connect(self._auto_bake)

        self._settings_agent_max_slope = QtWidgets.QDoubleSpinBox()
        self._settings_agent_max_slope.setValue(45.0)
        self._settings_agent_max_slope.setSingleStep(5.0)
        self._settings_agent_max_slope.setMinimum(0.1)
        self._settings_agent_max_slope.setMaximum(90.0)
        self._settings_agent_max_slope.valueChanged.connect(self._auto_bake)

        self._settings_region_min_size = QtWidgets.QSpinBox()
        self._settings_region_min_size.setValue(8)
        self._settings_region_min_size.setMinimum(1)
        self._settings_region_min_size.setMaximum(2048)
        self._settings_region_min_size.valueChanged.connect(self._auto_bake)

        self._settings_region_merge_size = QtWidgets.QSpinBox()
        self._settings_region_merge_size.setValue(20)
        self._settings_region_merge_size.setMinimum(1)
        self._settings_region_merge_size.setMaximum(64)
        self._settings_region_merge_size.valueChanged.connect(self._auto_bake)

        self._settings_edge_max_len = QtWidgets.QDoubleSpinBox()
        self._settings_edge_max_len.setValue(12.0)
        self._settings_edge_max_len.setSingleStep(1.0)
        self._settings_edge_max_len.setMinimum(0.001)
        self._settings_edge_max_len.valueChanged.connect(self._auto_bake)

        self._settings_edge_max_error = QtWidgets.QDoubleSpinBox()
        self._settings_edge_max_error.setValue(1.3)
        self._settings_edge_max_error.setSingleStep(0.1)
        self._settings_edge_max_error.setMinimum(0.001)
        self._settings_edge_max_error.valueChanged.connect(self._auto_bake)

        self._settings_verts_per_poly = QtWidgets.QSpinBox()
        self._settings_verts_per_poly.setValue(6)
        self._settings_verts_per_poly.setMinimum(3)
        self._settings_verts_per_poly.setMaximum(64)
        self._settings_verts_per_poly.valueChanged.connect(self._auto_bake)

        self._settings_detail_sample_distance = QtWidgets.QDoubleSpinBox()
        self._settings_detail_sample_distance.setValue(6.0)
        self._settings_detail_sample_distance.setSingleStep(1.0)
        self._settings_detail_sample_distance.setMinimum(0.001)
        self._settings_detail_sample_distance.valueChanged.connect(self._auto_bake)

        self._settings_detail_sample_maximum_error = QtWidgets.QDoubleSpinBox()
        self._settings_detail_sample_maximum_error.setValue(1.0)
        self._settings_detail_sample_maximum_error.setSingleStep(0.1)
        self._settings_detail_sample_maximum_error.setMinimum(0.001)
        self._settings_detail_sample_maximum_error.valueChanged.connect(self._auto_bake)

        self._options_auto_update = QtWidgets.QCheckBox()
        self._options_auto_update.setChecked(False)
        self._options_auto_update.stateChanged.connect(self._auto_bake)

        self._generate_button = QtWidgets.QPushButton("Bake navmesh")
        self._generate_button.clicked.connect(self._bake_navmesh)

        self._raster_layout = QtWidgets.QFormLayout()
        self._raster_layout.addRow("&Cell Size:", self._settings_cell_size)
        self._raster_layout.addRow("&Cell Height:", self._settings_cell_height)

        self._agent_layout = QtWidgets.QFormLayout()
        self._agent_layout.addRow("&Agent Height:", self._settings_agent_height)
        self._agent_layout.addRow("&Agent Radius:", self._settings_agent_radius)
        self._agent_layout.addRow("&Agent Max Climp:", self._settings_agent_max_climb)
        self._agent_layout.addRow("&Agent Max Slope:", self._settings_agent_max_slope)

        self._meshing_layout = QtWidgets.QFormLayout()
        self._meshing_layout.addRow("&Region Min Size:", self._settings_region_min_size)
        self._meshing_layout.addRow("&Region Merge Size:", self._settings_region_merge_size)
        self._meshing_layout.addRow("&Edge Max Length:", self._settings_edge_max_len)
        self._meshing_layout.addRow("&Edge Max Error:", self._settings_edge_max_error)
        self._meshing_layout.addRow("&Vertices Per Poly:", self._settings_verts_per_poly)
        self._meshing_layout.addRow("&Detail Sample Distance:", self._settings_detail_sample_distance)
        self._meshing_layout.addRow("&Detail Sample Max Error:", self._settings_detail_sample_maximum_error)

        raster_box = QtWidgets.QGroupBox("Rasterisation")
        agent_box = QtWidgets.QGroupBox("Agent")
        meshing_box = QtWidgets.QGroupBox("Meshing")

        raster_box.setLayout(self._raster_layout)
        agent_box.setLayout(self._agent_layout)
        meshing_box.setLayout(self._meshing_layout)

        self._layout = QtWidgets.QVBoxLayout()
        self._layout.addWidget(raster_box)
        self._layout.addWidget(agent_box)
        self._layout.addWidget(meshing_box)

        upd_layout = QtWidgets.QHBoxLayout()
        upd_layout.addWidget(self._options_auto_update)
        upd_label = QtWidgets.QLabel("Auto update")
        upd_layout.addWidget(upd_label, QtCore.Qt.AlignCenter)

        self._layout.addLayout(upd_layout)
        self._layout.addWidget(self._generate_button)

        self._layout.addStretch()

        self.setLayout(self._layout)

    def _auto_bake(self) -> None:
        if self._options_auto_update.isChecked():
            self._bake_navmesh()

    def _bake_navmesh(self) -> None:
        if self._navmesh is not None and self._navmesh_mode == 0:
            start_time: float = time.time()
            # get settings
            settings: Dict[str, Any] = self._navmesh.get_settings()
            settings["cellSize"] = self._settings_cell_size.value()
            settings["cellHeight"] = self._settings_cell_height.value()
            settings["agentHeight"] = self._settings_agent_height.value()
            settings["agentRadius"] = self._settings_agent_radius.value()
            settings["agentMaxClimb"] = self._settings_agent_max_climb.value()
            settings["agentMaxSlope"] = self._settings_agent_max_slope.value()
            settings["regionMinSize"] = self._settings_region_min_size.value()
            settings["regionMergeSize"] = self._settings_region_merge_size.value()
            settings["edgeMaxLen"] = self._settings_edge_max_len.value()
            settings["edgeMaxError"] = self._settings_edge_max_error.value()
            settings["vertsPerPoly"] = self._settings_verts_per_poly.value()
            settings["detailSampleDist"] = self._settings_detail_sample_distance.value()
            settings["detailSampleMaxError"] = self._settings_detail_sample_maximum_error.value()
            self._navmesh.set_settings(settings)
            self._navmesh.build_navmesh()

            finish_time: float = time.time()

            (vertices, polygons, sizes) = self._navmesh.get_navmesh_poligonization()
            self._canvas.draw_navmesh(vertices, polygons, sizes)

            self._main.set_status_message("Bake time: " + "{:.2f}".format(finish_time - start_time) + " sec. Vertices: " + str(len(vertices) // 3) + ". Polygons: " + str(len(sizes)))
        else:
            self._main.set_status_message("Before build you should load geometry from obj file")

    def load_level_command(self, obj_path: str) -> None:
        self._navmesh = rd.Navmesh()
        self._navmesh_mode = 0  # level
        self._navmesh.init_by_obj(obj_path)

    def load_bin_navmesh_command(self, bin_path: str) -> None:
        # at first clear canvas scene and create new navmesh instance
        self._navmesh = rd.Navmesh()
        self._navmesh_mode = 1  # binary data
        self._canvas.clear_scene()

        # next we should
        # 1. Init simple geometry
        # 2. Bake simple navmesh
        # 3. And only then load it instead baked mesh
        self._navmesh = rd.Navmesh()
        self._navmesh.init_by_raw([4.0, 0.0, 4.0, -4.0, 0.0, 4.0, -4.0, 0.0, -4.0, 4.0, 0.0, -4.0], [4, 0, 3, 2, 1])
        self._navmesh.build_navmesh()
        self._navmesh.load_navmesh(bin_path)
        (vertices, polygons, sizes) = self._navmesh.get_navmesh_poligonization()
        # draw navmesh into canvas and frame it, because the level is clear and camera move can be locked
        self._canvas.draw_navmesh(vertices, polygons, sizes, is_frame = True)
        self._main.set_status_message("Load navmesh from binary. Vertices: " + str(len(vertices) // 3) + ". Polygons: " + str(len(sizes)))

    def export_bin_navmesh_command(self, path: str) -> None:
        if self._navmesh is not None:
            self._navmesh.save_navmesh(path)
        else:
            self._main.set_status_message("Build navmesh before export")

    def export_txt_navmesh_command(self, path: str) -> None:
        if self._navmesh is not None:
            # write to the text file three strings: with vertex positions, with polygon corners and with polygon sizes
            # get polygonization
            (vertices, polygons, sizes) = self._navmesh.get_navmesh_poligonization()

            with open(path, "w") as file:
                file.write(" ".join([str(v) for v in vertices]))
                file.write("\n")
                file.write(" ".join([str(v) for v in polygons]))
                file.write("\n")
                file.write(" ".join([str(v) for v in sizes]))
        else:
            self._main.set_status_message("Build navmesh before export")


class ViewMainWidget(QtWidgets.QWidget):
    def __init__(self, main_window) -> None:
        super(ViewMainWidget, self).__init__()
        self._main_app = main_window
        self._layout = QtWidgets.QHBoxLayout(self)
        self._canvas_widget: CanvasWidget = CanvasWidget(main_window)
        self._settings_widget: SettingsWidget = SettingsWidget(self._canvas_widget, self)
        self._layout.addWidget(QtWidgets.QWidget.createWindowContainer(self._canvas_widget))
        self._layout.addWidget(self._settings_widget)

    def load_level(self, path: str) -> None:
        self._canvas_widget.load_level_command(path)
        self._settings_widget.load_level_command(path)

    def load_bin_navmesh(self, path: str) -> None:
        self._settings_widget.load_bin_navmesh_command(path)

    def export_bin_navmesh(self, path: str) -> None:
        self._settings_widget.export_bin_navmesh_command(path)

    def export_txt_navmesh(self, path: str) -> None:
        self._settings_widget.export_txt_navmesh_command(path)

    def set_status_message(self, message: str) -> None:
        self._main_app.set_status_message(message)


class ViewApp(QtWidgets.QMainWindow):
    def __init__(self) -> None:
        super(ViewApp, self).__init__()
        self._main: ViewMainWidget = ViewMainWidget(self)
        self.setCentralWidget(self._main)
        self.setWindowTitle("RecastDetour Viewer")
        self._status = self.statusBar()

        exit_action = QtGui.QAction("&Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.setStatusTip("Exit application")
        exit_action.triggered.connect(self.close)
        # load actions
        load_action = QtGui.QAction("&Load Level", self)
        load_action.setShortcut("Ctrl+L")
        load_action.setStatusTip("Load level obj file")
        load_action.triggered.connect(self._load_action_call)

        load_nm_action = QtGui.QAction("&Load Navmesh", self)
        load_nm_action.setShortcut("Ctrl+O")
        load_nm_action.setStatusTip("Load navmesh from binary file")
        load_nm_action.triggered.connect(self._load_nm_action_call)

        # export actions
        export_txt_action = QtGui.QAction("&Export txt...", self)
        export_txt_action.setShortcut("Ctrl+T")
        export_txt_action.setStatusTip("Export navmesh to txt file")
        export_txt_action.triggered.connect(self._export_txt_action_call)
        export_bin_action = QtGui.QAction("&Export bin...", self)
        export_bin_action.setShortcut("Ctrl+B")
        export_bin_action.setStatusTip("Export navmesh to binary file")
        export_bin_action.triggered.connect(self._export_bin_action_call)

        menubar = self.menuBar()
        file_menu = menubar.addMenu("&File")
        file_menu.addAction(load_action)
        file_menu.addAction(load_nm_action)

        export_submenu = file_menu.addMenu("&Export Navmeh")
        export_submenu.addAction(export_txt_action)
        export_submenu.addAction(export_bin_action)

        file_menu.addSeparator()
        file_menu.addAction(exit_action)

    def set_status_message(self, message: str) -> None:
        self._status.showMessage(message)

    def _load_action_call(self) -> None:
        dialog = QtWidgets.QFileDialog(self)
        dialog.setFileMode(QtWidgets.QFileDialog.ExistingFiles)
        dialog.setNameFilters(["Wavefront OBJ (*.obj)"])
        dialog.setViewMode(QtWidgets.QFileDialog.Detail)
        if dialog.exec():
            selected_files: List[str] = dialog.selectedFiles()
            if len(selected_files) > 0:
                selected_path: str = selected_files[0]
                ext: str = os.path.splitext(selected_path)[1]
                if ext == ".obj":
                    self._main.load_level(selected_path)

    def _load_nm_action_call(self) -> None:
        dialog = QtWidgets.QFileDialog(self)
        dialog.setFileMode(QtWidgets.QFileDialog.ExistingFiles)
        dialog.setNameFilters(["RecastDetour Binary (*.bin)"])
        dialog.setViewMode(QtWidgets.QFileDialog.Detail)
        if dialog.exec():
            selected_files: List[str] = dialog.selectedFiles()
            if len(selected_files) > 0:
                selected_path: str = selected_files[0]
                ext: str = os.path.splitext(selected_path)[1]
                if ext == ".bin":
                    self._main.load_bin_navmesh(selected_path)

    def _export_txt_action_call(self) -> None:
        dialog_result = QtWidgets.QFileDialog.getSaveFileName(self, "Export Text", "", "Text File (*.txt)")
        file_path: str = dialog_result[0]
        if len(dialog_result) > 0:
            file_path = os.path.splitext(file_path)[0] + ".txt"
            self._main.export_txt_navmesh(file_path)

    def _export_bin_action_call(self) -> None:
        dialog_result = QtWidgets.QFileDialog.getSaveFileName(self, "Export Binary", "", "RecastDetour Binary (*.bin)")
        file_path: str = dialog_result[0]
        if len(dialog_result) > 0:
            file_path = os.path.splitext(file_path)[0] + ".bin"
            self._main.export_bin_navmesh(file_path)


if __name__ == "__main__":
    app = QtWidgets.QApplication()
    main: ViewApp = ViewApp()
    main.resize(1200, 820)
    main.show()
    app.exec()
