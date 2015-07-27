#include "sphere.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/**********************************************************************
 * This function intersects a ray with a given sphere 'sph'. You should
 * use the parametric representation of a line and do the intersection.
 * The function should return the parameter value for the intersection, 
 * which will be compared with others to determine which intersection
 * is closest. The value -1.0 is returned if there is no intersection
 *
 * If there is an intersection, the point of intersection should be
 * stored in the "hit" variable
 **********************************************************************/

// Function to calculate the components for caluclating t
// Defined as a separate function so it can be used for determining shadows as well
Vector intersect_sphere_components(Point o, Vector u, Spheres *sph){

  // From "12. Ray Tracing.pdf" slide # 33
  float A, B, C, radiusPow2;
  Vector sphCentre;

  // A = Xd^2 + Yd^2 + Zd^2
  A = pow(u.x, 2) + pow(u.y, 2) + pow(u.z, 2);

  // Define X0, Y0, Z0
  sphCentre.x = sph->center.x;
  sphCentre.y = sph->center.y;
  sphCentre.z = sph->center.z;

  // B = 2(Xd(Xs - X0) + Yd(Ys - Y0) + Zd(Zs-Z0))
  B = 2*( u.x*(o.x - sphCentre.x) + u.y*(o.y - sphCentre.y) + u.z*(o.z - sphCentre.z) );

  // Define r^2
  radiusPow2 = pow(sph->radius, 2);

  // C = (Xs - X0)^2 + (Ys-Y0)^2 + (Zs - Z0)^ 2 - r^2
  C = pow((o.x - sphCentre.x), 2) + pow((o.y - sphCentre.y), 2) + pow((o.z - sphCentre.z), 2) - radiusPow2;

  Vector components = {A, B, C};

  return components;

}


float intersect_sphere(Point o, Vector u, Spheres *sph, Point *hit) {

  Vector components = intersect_sphere_components(o, u, sph);

  float t, t1, t2, A, B, C, squareRoot;

  A = components.x;
  B = components.y;
  C = components.z;
  squareRoot = (pow(B, 2) - 4*A*C);

  // t(1,2) = (-B (+,-) sqrt(B^2 - 4AC))/2A
  // We need to make sure the sqrt() part is not negative in order to avoid errors
  if (squareRoot < 0){
    return -1.0;
  }

  else{
    // Calculate t1 and t2 now
    t1 = (-B + sqrt(squareRoot))/(2*A);
    t2 = (-B - sqrt(squareRoot))/(2*A);

    // Use the smallest value since it's closer
    t = t1;
    if (t1 > t2 && t2 > 0.02){ // Check with a value greater than 0 to avoid self-shadows
      t = t2;
    }
  
    // Make sure it is not negative
    if (t < 0.02){ // Check with a value greater than 0 to avoid self-shadows
      return -1.0;
    }

    // If there is an intersection, the point of intersection should be
    // stored in the "hit" variable
    // From "12. Ray Tracing.pdf" slide # 23
    hit->x = o.x + t2*u.x;
    hit->y = o.y + t2*u.y;
    hit->z = o.z + t2*u.z;

    return t;
  }

}

/*********************************************************************
 * This function returns a pointer to the sphere object that the
 * ray intersects first; NULL if no intersection. You should decide
 * which arguments to use for the function. For exmaple, note that you
 * should return the point of intersection to the calling function.
 **********************************************************************/
// Using the same parameters as the previous function
Spheres *intersect_scene(Point o, Vector u, Spheres *sph, Point *hit) {

  // Set up variables
  Spheres *sphereInput = sph;
  Spheres *nearestSphere = NULL;// NULL since we have not detected an intersection yet
  float distance, intersectionPoint;
  bool isCloser, PosIntersect;
  distance = FLT_MAX;   // The farthest distance initially

  while (sphereInput != NULL){

    // Calculate intersection
    intersectionPoint = intersect_sphere(o, u, sphereInput, hit);

    // Intersection conditions
    isCloser = distance > intersectionPoint;
    PosIntersect = (intersectionPoint != -1.0);

    // Check if it intersects
    if ( isCloser && PosIntersect ){
      distance = intersectionPoint;
      nearestSphere = sphereInput;
    }
    // Move on to the next sphere
    sphereInput = sphereInput->next;
  }
  return nearestSphere;
}

/*****************************************************
 * This function adds a sphere into the sphere list
 *
 * You need not change this.
 *****************************************************/
Spheres *add_sphere(Spheres *slist, Point ctr, float rad, float amb[],
		    float dif[], float spe[], float shine, 
		    float refl, int sindex) {
  Spheres *new_sphere;

  new_sphere = (Spheres *)malloc(sizeof(Spheres));
  new_sphere->index = sindex;
  new_sphere->center = ctr;
  new_sphere->radius = rad;
  (new_sphere->mat_ambient)[0] = amb[0];
  (new_sphere->mat_ambient)[1] = amb[1];
  (new_sphere->mat_ambient)[2] = amb[2];
  (new_sphere->mat_diffuse)[0] = dif[0];
  (new_sphere->mat_diffuse)[1] = dif[1];
  (new_sphere->mat_diffuse)[2] = dif[2];
  (new_sphere->mat_specular)[0] = spe[0];
  (new_sphere->mat_specular)[1] = spe[1];
  (new_sphere->mat_specular)[2] = spe[2];
  new_sphere->mat_shineness = shine;
  new_sphere->reflectance = refl;
  new_sphere->next = NULL;

  if (slist == NULL) { // first object
    slist = new_sphere;
  } else { // insert at the beginning
    new_sphere->next = slist;
    slist = new_sphere;
  }

  return slist;
}

/******************************************
 * computes a sphere normal - done for you
 ******************************************/
Vector sphere_normal(Point q, Spheres *sph) {
  Vector rc;

  rc = get_vec(sph->center, q);
  normalize(&rc);
  return rc;
}


// Function to see if the current sphere is in a shadow or not
bool intersect_sphere_shadow(Point q, Vector u, Spheres *sph){

  while (sph != NULL){

    
    float t1, t2, A, B, C, squareRoot;

    // compute A, B, C
    normalize(&u);
    Vector components = intersect_sphere_components(q, u, sph);
    A = components.x;
    B = components.y;
    C = components.z;

    // Compute t1 and t2
    squareRoot = (pow(B, 2) - 4*A*C);
    t1 = (-B + sqrt(squareRoot))/(2*A);
    t2 = (-B - sqrt(squareRoot))/(2*A);

    // Check if the squre root part o
    if (squareRoot > 0){
      if ((t1 > 0.02) || (t2 > 0.02)) // Check with a value greater than 0 to avoid self-shadows
        return true;
    }
    sph = sph->next;

  }
  return false;
}