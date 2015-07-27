#include <stdio.h>
#include <GL/glut.h>
#include <math.h>
#include "global.h"
#include "sphere.h"

//
// Global variables
//
extern int win_width;
extern int win_height;

extern GLfloat frame[WIN_HEIGHT][WIN_WIDTH][3];  

extern float image_width;
extern float image_height;

extern Point eye_pos;
extern float image_plane;
extern RGB_float background_clr;
extern RGB_float null_clr;

extern Spheres *scene;

// light 1 position and color
extern Point light1;
extern float light1_ambient[3];
extern float light1_diffuse[3];
extern float light1_specular[3];

// global ambient term
extern float global_ambient[3];

// light decay parameters
extern float decay_a;
extern float decay_b;
extern float decay_c;

extern int shadow_on;
extern int step_max;
extern int reflect_on;
extern int chess_on;
extern int stochast_on;

/////////////////////////////////////////////////////////////////////

/*********************************************************************
 * Phong illumination - you need to implement this!
 *********************************************************************/
RGB_float phong(Point q, Vector v, Vector surf_norm, Spheres *sph) {

  RGB_float color = {0, 0, 0}; //initialize

  // Phong's local illumination model
  // Global ambient + ambient + decay*(diffiuse + specular)
  // I = Iga * Kga + Ia * Ka + (1/(a+bd+cd^2)) * (Id * Kd *(n*l) + Is * Ks *(r*v)^N)
  // I = Iga * Kga + Ia * Ka + (1/(a+bd+cd^2))*(Id * Kd *(n*l)) + (1/(a+bd+cd^2))*(Is * Ks *(r*v)^N)

  // First turn the two points q and the light source into a vector
  Vector lightRay = get_vec(q, light1);
  normalize(&lightRay);

  // "d is the distance between the light source and the point on the object"
  float d = vec_len(lightRay);


  // Global ambient: Iga * Kga
  // Apply thevalues contained in the global ambient array to the sphere
  /*color.r = color.r + global_ambient[0]*sph->reflectance;
  color.g = color.g + global_ambient[1]*sph->reflectance;
  color.b = color.b + global_ambient[2]*sph->reflectance;*/
  RGB_float ambient = {0, 0, 0};
  ambient.r = ambient.r + global_ambient[0]*sph->reflectance;
  ambient.g = ambient.g + global_ambient[1]*sph->reflectance;
  ambient.b = ambient.b + global_ambient[2]*sph->reflectance;


  // Ambient: Ia * Ka
  ambient.r = ambient.r + light1_ambient[0]*sph->mat_ambient[0];
  ambient.g = ambient.g + light1_ambient[1]*sph->mat_ambient[1];
  ambient.b = ambient.b + light1_ambient[2]*sph->mat_ambient[2];


  // 1/(a+bd+cd^2)
  float light_decay = 1 / (decay_a + decay_b*d + decay_c*pow(d,2) );


  // decay*Diffuse: (1/(a+bd+cd^2))*(Id * Kd *(n*l))
  float lightNormVal = vec_dot(surf_norm, lightRay); // (n*l)
  color.r = color.r + light_decay*(light1_diffuse[0] * sph->mat_diffuse[0] * lightNormVal);
  color.g = color.g + light_decay*(light1_diffuse[1] * sph->mat_diffuse[1] * lightNormVal);
  color.b = color.b + light_decay*(light1_diffuse[2] * sph->mat_diffuse[2] * lightNormVal);


  // decay*Specular: (1/(a+bd+cd^2))*(Is * Ks *(r*v)^N)
  // compute r first
  float angle = vec_dot(surf_norm, lightRay);
  if (angle < 0) angle = 0;

  Vector scaledNorm = vec_scale(surf_norm, 2*angle);
  Vector reflec = vec_minus(scaledNorm, lightRay);
  normalize(&reflec);

  // compute (r*v)^N
  float reflecViewVal = vec_dot(reflec, v);
  reflecViewVal = pow(reflecViewVal, sph->mat_shineness);

  // compute decay*specular
  // lightDecay*(Is * Ks * (r*v)^N)
  color.r = color.r + light_decay*(light1_specular[0] * sph->mat_specular[0] * reflecViewVal);
  color.g = color.g + light_decay*(light1_specular[1] * sph->mat_specular[1] * reflecViewVal);
  color.b = color.b + light_decay*(light1_specular[2] * sph->mat_specular[2] * reflecViewVal);

  // Check if shadows are enabled
  bool intersectShadow = intersect_sphere_shadow(q, lightRay, scene);
  if (shadow_on && intersectShadow){
    color = ambient;
  }

	return color;
}

/////////////////////////REFLECTION/////////////////////////////////

// Function for computing the color of the reflected color
Vector reflect_vector(Vector v, Vector norm){

  normalize(&v);
  float scale = -2 * vec_dot(norm, v);
  Vector scaledNorm = vec_scale(norm, scale);

  // Subtract the scaled surface norm and the light vector
  Vector vec = vec_plus(scaledNorm, v);
  normalize(&vec);

  return vec;

}

/////////////////////CHESSBOARD/////////////////////////////////////

// From "12. Ray Tracing.pdf" slide # 34
Vector chess_surf_norm = {0, 1, 1};
Point chess_position = {1, 1, 1};

int chessboard(Point q, Vector v, Point *hit){

  // Calculate t = -(AXs + BYs + CZs + D)/(AXd +BYd + CZd)
  float fraction1, fraction2, t;

  normalize(&chess_surf_norm);

  //-(AXs + BYs + CZs + D)
  fraction1 = -1.0 * vec_dot(chess_surf_norm, get_vec(q,chess_position));

  //(AXd +BYd + CZd)
  fraction2 = vec_dot(chess_surf_norm, v);

  // t = -(AXs + BYs + CZs + D)/(AXd +BYd + CZd)
  t = (fraction1/fraction2);

  // if t = 0 then there is not an intersection
  if (fraction2 != 0 && t==0){
    return 0;
  }

  // if t > 0 then there is an intersection
  else if (t > 0){
    hit->x = q.x + (v.x * t);
    hit->y = q.y + (v.y * t);
    hit->z = q.z + (v.z * t);
    return 1;
  }

  else
    return 0;
}

RGB_float color_chessboard(Point *hit){

// round to use modular to estimate positions of hit 
  int roundedChessHitx = int(hit->x);
  int roundedChessHitz = int(hit->z);

  RGB_float color;


  if((roundedChessHitx % 2) == 0){

    if ((roundedChessHitz % 2) == 0){
      color.r = 0;
      color.g = 0;
      color.b = 0;
    }

    else if ((roundedChessHitz % 2) != 0){
      color.r = 1;
      color.g = 1;
      color.b = 1;
    }  
  }

  else if ((roundedChessHitx % 2) != 0){
    if ((roundedChessHitz % 2) == 0){
      color.r = 1;
      color.g = 1;
      color.b = 1;
    }
  }

  // Limit it to 8x8
  if ((roundedChessHitx > 4) || (roundedChessHitx <-4) ||
      (roundedChessHitz > 4) || (roundedChessHitz <-4))
      
      color = background_clr;

  return color;
}


/************************************************************************
 * This is the recursive ray tracer - you need to implement this!
 * You should decide what arguments to use.
 ************************************************************************/
RGB_float recursive_ray_trace(int stepIndex, Point q, Vector v) {

	// Initialize the colour to background in case it doesn't detect any objects
  RGB_float color = background_clr;
  RGB_float reflec_color = {0, 0, 0,}; // Initialize

  // Use the intersect scene function from sphere.cpp
  Point *hit = new Point;
  Point *chessHit = new Point;
  Spheres *nearestSph;
  nearestSph = intersect_scene(q, v, scene, hit);
    
  // Chessboard computation
  if (chess_on && (chessboard(q, v, chessHit) == 1)){
      
    Vector eye = get_vec(*chessHit, eye_pos);
    normalize(&eye);

    color = color_chessboard(chessHit);

    // Chessboard shadow
    Vector shade = get_vec(*chessHit, light1);
    Spheres *nullSphere = NULL;

    if ((intersect_sphere_shadow(*chessHit, shade, nullSphere)) && shadow_on)
      color = clr_scale(color, 0.5);

  }

  // Check if the light rays hit anything
  if (nearestSph != NULL){

    // Configure the inputs for phong()
    Vector eye = get_vec(*hit, eye_pos);
    normalize(&eye);

    Vector surf_norm = sphere_normal(*hit, nearestSph);

    color = phong(*hit, eye, surf_norm, nearestSph);


    // Relfection computation
    // Check if reflection is enabled and if the number of current reflections are less than the maximum
    if (reflect_on && (step_max > stepIndex)){

      Vector reflection_vec = reflect_vector(v, surf_norm);
      stepIndex = stepIndex + 1;

      // Use recursion
      reflec_color = recursive_ray_trace(stepIndex, *hit, reflection_vec);

      // Stochastic ray gen
      if (stochast_on == 1){

        // Duplicate reflection_vec to modify it
        Vector stochastv = reflection_vec;

        // Generate random rays
        for (int i = 0; i < STOCHASTIC_RAYS_COUNT; i++){

          stochastv.x = stochastv.x + rand();
          stochastv.y = stochastv.y + rand();
          stochastv.z = stochastv.z + rand();

          normalize(&stochastv);

          // Apply the randomly generated values to the color output
          RGB_float randomColor = recursive_ray_trace(stepIndex, *hit, stochastv);
          randomColor = clr_scale(randomColor, 0.02);
          reflec_color = clr_add(reflec_color, randomColor);
        }
      }
      
      reflec_color = clr_scale(reflec_color, nearestSph->reflectance);
      color = clr_add(color, reflec_color);
    }

  }

	return color;
}

/*********************************************************************
 * This function traverses all the pixels and cast rays. It calls the
 * recursive ray tracer and assign return color to frame
 *
 * You should not need to change it except for the call to the recursive
 * ray tracer. Feel free to change other parts of the function however,
 * if you must.
 *********************************************************************/
void ray_trace() {
  int i, j;
  float x_grid_size = image_width / float(win_width);
  float y_grid_size = image_height / float(win_height);
  float x_start = -0.5 * image_width;
  float y_start = -0.5 * image_height;
  RGB_float ret_color;
  Point cur_pixel_pos;
  Vector ray;

  // ray is cast through center of pixel
  cur_pixel_pos.x = x_start + 0.5 * x_grid_size;
  cur_pixel_pos.y = y_start + 0.5 * y_grid_size;
  cur_pixel_pos.z = image_plane;

  for (i=0; i<win_height; i++) {
    for (j=0; j<win_width; j++) {
      ray = get_vec(eye_pos, cur_pixel_pos);

      //
      // You need to change this!!!
      //
      ret_color = recursive_ray_trace(0, eye_pos, ray);

      // ret_color = background_clr; // just background for now

      // Parallel rays can be cast instead using below
      //
      // ray.x = ray.y = 0;
      // ray.z = -1.0;
      // ret_color = recursive_ray_trace(cur_pixel_pos, ray, 1);

      // Checkboard for testing
      RGB_float clr = {float(i/32), 0, float(j/32)};
     // ret_color = clr;

      frame[i][j][0] = GLfloat(ret_color.r);
      frame[i][j][1] = GLfloat(ret_color.g);
      frame[i][j][2] = GLfloat(ret_color.b);

      cur_pixel_pos.x += x_grid_size;
    }

    cur_pixel_pos.y += y_grid_size;
    cur_pixel_pos.x = x_start;
  }
}
