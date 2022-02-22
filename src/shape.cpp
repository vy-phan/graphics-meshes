#include "shape.h"

#include "types.h"

#include "debugging.h"

#include <cmath> // std::sin(), std::cos()

using namespace glm;

namespace graphics101 {

void Cylinder::tessellate( Mesh& mesh )
{
    mesh.clear();

    real r = 1.; // radius = 1

    // Iterate around the cylinder
    for ( int i = 0; i < m_slices; ++i ) {
        // u = i / slices
        // phi = 2*pi*u
        real u = real(i) / real( m_slices ) ;
        real phi = 2. * pi * u;

        // x = r*cos(phi)
        // y = r*sin(phi)
        real x = r*glm::cos( phi );
        real y = r*glm::sin( phi );

        // Add the normal of the position at x,y
        mesh.normals.push_back( normalize( vec3( x, y, 0 ) ) );

        for ( int j = 0; j < m_stacks; ++j ) {
            real z = real(j) / real( m_stacks - 1 );

            // Add x,y,z Cartesian coordinates to mesh's positions
            mesh.positions.push_back( vec3( x, y, z ) );

            // Add u,v texture coordinates to mesh's texcord
            mesh.texcoords.push_back( vec2( u, z ) );

            if ( i > 0 && j > 0 ) {
                // Triangle 1 - [ upper-left, lower-left, lower-right ]
                mesh.face_positions.push_back( Triangle( i*m_stacks + j,
                                                         (i-1)*m_stacks + j,
                                                         (i-1)*m_stacks + (j - 1) ) );
                mesh.face_texcoords.push_back( Triangle( i*m_stacks + j,
                                                         (i-1)*m_stacks + j,
                                                         (i-1)*m_stacks + (j - 1) ) );
                mesh.face_normals.push_back( Triangle( i, i-1, i-1 ) );

                // Triangle 2 - [ upper-left, upper-right, lower-right ]
                mesh.face_positions.push_back( Triangle( i*m_stacks + j,
                                                         (i-1)*m_stacks + (j - 1),
                                                         i*m_stacks + (j - 1) ) );
                mesh.face_texcoords.push_back( Triangle( i*m_stacks + j,
                                                         (i-1)*m_stacks + (j - 1),
                                                         i*m_stacks + (j - 1) ) );
                mesh.face_normals.push_back( Triangle( i, i-1, i ) );
            }

        }

    }

    //Stick the last edge
    for ( int j = 0; j < m_stacks; ++j ) {
        // Add extra 1,v texture coordinates to mesh's texcord
        mesh.texcoords.push_back( vec2( 1.0, real(j) / real( m_stacks - 1 ) ) );

        if ( j > 0 ) {
            // Triangle 1 - [ upper-right, upper-left, lower-left ]
            mesh.face_positions.push_back( Triangle( j,
                                                     (m_slices-1)*m_stacks + j,
                                                     (m_slices-1)*m_stacks + (j - 1) ) );
            mesh.face_texcoords.push_back( Triangle( m_slices*m_stacks + j,
                                                     (m_slices-1)*m_stacks + j,
                                                     (m_slices-1)*m_stacks + (j - 1) ) );
            mesh.face_normals.push_back( Triangle( 0, m_slices-1, m_slices-1 ) );

            // Triangle 2 - [ upper-right, lower-left, lower-right ]
            mesh.face_positions.push_back( Triangle( j,
                                                     (m_slices-1)*m_stacks + (j - 1),
                                                     j-1 ) );
            mesh.face_texcoords.push_back( Triangle( m_slices*m_stacks + j,
                                                     (m_slices-1)*m_stacks + (j - 1),
                                                     m_slices*m_stacks + (j - 1) ) );
            mesh.face_normals.push_back( Triangle( 0, m_slices-1, 0 ) );
        }
    }

    // Assert
    assert( mesh.positions.size() == m_slices*m_stacks );
    assert( mesh.texcoords.size() == ( m_slices + 1 )*m_stacks );
    assert( mesh.normals.size() == m_slices );

}

void Sphere::tessellate( Mesh& mesh )
{
    mesh.clear();

    real r = 1.; // radius = 1

    for ( int i = 0; i < m_slices; ++i ) {
        // u = i / slices
        // phi = 2*pi*u
        real u = ( real(i) / real( m_slices ) );
        real phi = 2. * pi * u;


        for ( int j = 0; j < m_stacks; ++j ) {
            real v = real(j + 1) / real( m_stacks + 2 - 1 );
            real theta = pi * v;

            real x = r*glm::sin( theta )*glm::cos( phi );
            real y = r*glm::sin( theta )*glm::sin( phi );
            real z = r*glm::cos( -1 * theta );

            // Add x,y,z Cartesian coordinates to mesh's positions
            mesh.positions.push_back( vec3( x, y, z ) );

            // Add u,v texture coordinates to mesh's texcord
            mesh.texcoords.push_back( vec2( u, 1 - v ) );

            // Add the normal of the position at x,y
            mesh.normals.push_back( normalize( vec3( x, y, z ) ) );

            if ( i > 0 && j > 0 ) {
                // Triangle 1 - [ upper-left, lower-left, lower-right ]
                mesh.face_positions.push_back( Triangle( i*m_stacks + j,
                                                         (i-1)*m_stacks + (j - 1),
                                                         (i-1)*m_stacks + j ) );
                mesh.face_texcoords.push_back( Triangle( i*m_stacks + j,
                                                         (i-1)*m_stacks + (j - 1),
                                                         (i-1)*m_stacks + j ) );
                mesh.face_normals.push_back( Triangle( i*m_stacks + j,
                                                       (i-1)*m_stacks + (j - 1),
                                                       (i-1)*m_stacks + j ) ) ;

                // Triangle 2 - [ upper-left, upper-right, lower-right ]
                mesh.face_positions.push_back( Triangle( i*m_stacks + j,
                                                         i*m_stacks + (j - 1),
                                                         (i-1)*m_stacks + (j - 1) ) );
                mesh.face_texcoords.push_back( Triangle( i*m_stacks + j,
                                                         i*m_stacks + (j - 1),
                                                         (i-1)*m_stacks + (j - 1) ) );
                mesh.face_normals.push_back( Triangle( i*m_stacks + j,
                                                       i*m_stacks + (j - 1),
                                                       (i-1)*m_stacks + (j - 1) )  );
            }

        }

    }

    //Stick the last edge
    for ( int j = 0; j < m_stacks; ++j ) {
        // Add extra 1,v texture coordinates to mesh's texcord
        mesh.texcoords.push_back( vec2( 1.0, 1 - ( real(j+1) / real( m_stacks + 2 - 1 ) ) ) );

        if ( j > 0 ) {
            // Triangle 1 - [ upper-right, upper-left, lower-left ]
            mesh.face_positions.push_back( Triangle( j,
                                                     (m_slices-1)*m_stacks + (j - 1),
                                                     (m_slices-1)*m_stacks + j ) );
            mesh.face_texcoords.push_back( Triangle( m_slices*m_stacks + j,
                                                     (m_slices-1)*m_stacks + (j - 1),
                                                     (m_slices-1)*m_stacks + j ) );
            mesh.face_normals.push_back( Triangle( j,
                                                   (m_slices-1)*m_stacks + (j - 1),
                                                   (m_slices-1)*m_stacks + j ) );

            // Triangle 2 - [ upper-right, lower-left, lower-right ]
            mesh.face_positions.push_back( Triangle( j,
                                                     j-1,
                                                     (m_slices-1)*m_stacks + (j - 1) ) );
            mesh.face_texcoords.push_back( Triangle( m_slices*m_stacks + j,
                                                     m_slices*m_stacks + (j - 1),
                                                     (m_slices-1)*m_stacks + (j - 1) ) );
            mesh.face_normals.push_back( Triangle( j,
                                                   j-1,
                                                   (m_slices-1)*m_stacks + (j - 1) ) );
        }
    }

    // Add poles
    mesh.positions.push_back( vec3( 0., 0., 1. ) ); // North Pole
    mesh.positions.push_back( vec3( 0., 0., -1. ) ); // South Pole

    mesh.texcoords.push_back( vec2( 0.5, 0. ) );
    mesh.texcoords.push_back( vec2( 0.5, 1. ) );

    mesh.normals.push_back( vec3( 0., 0., 1. ) );
    mesh.normals.push_back( vec3( 0., 0., -1. ) );

    // Index into poles
    int north_p = mesh.positions.size() - 2;
    int south_p = mesh.positions.size() - 1;
    int south_uv = mesh.texcoords.size() - 2;
    int north_uv = mesh.texcoords.size() - 1;

    // Stitch poles
    for ( int i = 0; i < m_slices; ++i ) {
        if ( i > 0 ) {
            // South pole
            mesh.face_positions.push_back( Triangle( south_p,
                                                     (i+1)*m_stacks - 1,
                                                     i*m_stacks - 1 ) );
            mesh.face_texcoords.push_back( Triangle( south_uv,
                                                     i*m_stacks - 1,
                                                     (i+1)*m_stacks - 1 ) );
            mesh.face_normals.push_back( Triangle( south_p,
                                                   (i+1)*m_stacks - 1,
                                                   i*m_stacks - 1 ) );
            // North pole
            mesh.face_positions.push_back( Triangle( north_p,
                                                     (i-1)*m_stacks,
                                                     i*m_stacks ) );
            mesh.face_texcoords.push_back( Triangle( north_uv,
                                                     (i-1)*m_stacks,
                                                     (i+1)*m_stacks ) );
            mesh.face_normals.push_back( Triangle( north_p,
                                                   (i-1)*m_stacks,
                                                   i*m_stacks ) );
        }
        else {
            // South pole
            mesh.face_positions.push_back( Triangle( south_p,
                                                     m_stacks-1,
                                                     (m_slices-1)*m_stacks - 2 ) );
            mesh.face_texcoords.push_back( Triangle( south_uv,
                                                     (m_slices-1)*m_stacks + m_stacks - 1,
                                                     m_slices*m_stacks + m_stacks - 1) );
            mesh.face_normals.push_back( Triangle( south_p,
                                                   m_stacks-1,
                                                   (m_slices-1)*m_stacks - 2 ) );
            // North pole
            mesh.face_positions.push_back( Triangle( north_p,
                                                     (m_slices-1)*m_stacks,
                                                     0 ) );
            mesh.face_texcoords.push_back( Triangle( north_uv,
                                                     (m_slices-1)*m_stacks,
                                                     m_slices*m_stacks ) );
            mesh.face_normals.push_back( Triangle( north_p,
                                                   (m_slices-1)*m_stacks,
                                                   0 ) );
        }
    }

    // Assert
    assert( mesh.positions.size() == m_slices*m_stacks + 2 );
    assert( mesh.texcoords.size() == ( m_slices + 1 )*m_stacks + 2);
    assert( mesh.normals.size() == m_slices*m_stacks + 2);
}

void Cone::tessellate( Mesh& mesh )
{
    mesh.clear();

    // Your code goes here.
}

void Torus::tessellate( Mesh& mesh )
{
    mesh.clear();

    // Your code goes here.
}

void Cube::tessellate( Mesh& mesh )
{
    mesh.clear();

    // Your code goes here.
}

}
