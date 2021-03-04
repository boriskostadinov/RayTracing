//#include "color.h"
//#include "ray.h"
//#include "vec3.h"


#include<fstream>
#include<cstdlib>
using namespace std;

#include "rtweekend.h"
#include "camera.h"
#include "color.h"
#include "hittable_list.h"
#include "material.h"
#include "moving_sphere.h"
#include "sphere.h"
#include "aarect.h"
#include "box.h"


#include <iostream>
#include <map>



#define NOMINMAX

// optix
#include "primeCommon.h"
#include <optixu/optixu_math_namespace.h>
#include <sutil.h>

/*
// TODO*(boris): Add comment and Math here
double hit_sphere(const point3& center, double radius, const ray& r) {
    vec3 oc = r.origin() - center;
    auto a = r.direction().length_squared();
    auto half_b = dot(oc, r.direction());
    auto c = oc.length_squared() - radius * radius;
    auto discriminant = half_b * half_b - a * c;

    if (discriminant < 0) {
        return -1.0;
    }
    else {
        return (-half_b - sqrt(discriminant)) / a;
    }
}
*/

auto material_global = make_shared<lambertian>(color(0.8, 0.8, 0.0));


color ray_color(const ray& r, const color& background, const hittable& world, int depth) {
    hit_record rec;

    // If we've exceeded the ray bounce limit, no more light is gathered.
    if (depth <= 0)
        return color(0, 0, 0);

/*  if (world.hit(r, 0.001, infinity, rec)) {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth - 1);
        return color(0, 0, 0);
    }

    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5 * (unit_direction.y() + 1.0);
    return (1.0 - t) * color(1.0, 1.0, 1.0) + t * color(0.5, 0.7, 1.0);
*/
    // If the ray hits nothing, return the background color.
    if (!world.hit(r, 0.001, infinity, rec))
        return background;

    ray scattered;
    color attenuation;
    color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

    if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered))
        return emitted;

    return emitted + attenuation * ray_color(scattered, background, world, depth - 1);

}

// vec3 -> float3 ; to OptiX
float3 make_float3(const vec3& vec) {
    float3 f;
    f.x = vec.x();
    f.y = vec.y();
    f.z = vec.z();
    return f;
}

// float3 -> vec3 ; from OptiX
vec3 make_vec3(const float3& fl3) {
    vec3 v;
    v.e[0] = fl3.x;
    v.e[1] = fl3.y;
    v.e[2] = fl3.z;

    return v;
}

// R to r
ray make_ray(const Ray& R) {
    ray r;
    r.origin() = make_vec3(R.origin);
    r.direction() = make_vec3(R.dir);

    return r;
}


color rt_color(Hit h, Ray r, const color& background, PrimeMesh mesh) {
    hit_record rec;
    rec.t = h.t;
    rec.u = h.u;
    rec.v = h.v;
    
    //float3 checker = make_float3(rec.p);
    //checker = r.origin + r.dir * h.t;
    //rec.p = r.origin + r.dir * h.t;  // NEVER CHANGE

    make_float3(rec.p) = r.origin + r.dir * h.t;

    vec3 re_checked = rec.p;

    int3* indices = mesh.getVertexIndices();
    float3* vertices = mesh.getVertexData();


    // calc normal and matherial
    if (h.t < 0.0f)
    {
        return background;
    }
    else
    {
        int3 tri = indices[h.triId];
        float3 v0 = vertices[tri.x];
        float3 v1 = vertices[tri.y];
        float3 v2 = vertices[tri.z];
        float3 e0 = v1 - v0;
        float3 e1 = v2 - v0;
        float3 n = optix::normalize(optix::cross(e0, e1));

        //rec.normal = n;
        rec.normal = make_vec3(n);
        rec.mat_ptr = material_global;

        //R->r TODO

        //ray ray = r;
        ray r1 = make_ray(r);

        ray scattered;
        color attenuation;
        color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

        if (!rec.mat_ptr->scatter(r1, rec, attenuation, scattered))
            return emitted;

        return emitted + attenuation;
    }

}





// A very complex scene
hittable_list random_scene() {
    hittable_list world;

    auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));
    world.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(checker)));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

            if ((center - point3(4, 0.2, 0)).length() > 0.9) {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    auto center2 = center + vec3(0, random_double(0, .5), 0);
                    world.add(make_shared<moving_sphere>(
                        center, center2, 0.0, 1.0, 0.2, sphere_material));
                }
                else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = make_shared<dielectric>(1.5);
    world.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    world.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));

    auto material3 = make_shared<metal>(color(0.7, 0.6, 0.5), 0.0);
    world.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    return world;
}

hittable_list two_spheres() {
    hittable_list objects;

    auto checker = make_shared<checker_texture>(color(0.2, 0.3, 0.1), color(0.9, 0.9, 0.9));

    objects.add(make_shared<sphere>(point3(0, -10, 0), 10, make_shared<lambertian>(checker)));
    objects.add(make_shared<sphere>(point3(0, 10, 0), 10, make_shared<lambertian>(checker)));

    return objects;
}

hittable_list earth() {
    auto earth_texture = make_shared<image_texture>("earthmap.jpg");
    auto earth_surface = make_shared<lambertian>(earth_texture);
    auto globe = make_shared<sphere>(point3(0, 0, 0), 2, earth_surface);

    return hittable_list(globe);
}

hittable_list simple_light() {
    hittable_list objects;

    auto pertext = make_shared<image_texture>("earthmap.jpg");
    objects.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(pertext)));
    objects.add(make_shared<sphere>(point3(0, 2, 0), 2, make_shared<lambertian>(pertext)));

    auto difflight = make_shared<diffuse_light>(color(8, 8, 0.5));
    //objects.add(make_shared<xy_rect>(3, 5, 1, 3, -2, difflight));
    objects.add(make_shared<sphere>(point3(0, 7, 0), 1, difflight));
    return objects;
}


int main() {
    // MeshFile
    std::string objFilename;
    bool meshIncluded = 0;

    // Image
    double aspect_ratio;
    int image_width;       
    int samples_per_pixel;    
    int max_depth;   
    color background;

    
    // Camera
    point3 lookfrom;    
    point3 lookat;    
    vec3 vup;
    double vfov;
    double aperture;
    double dist_to_focus;

    // Objects
    hittable_list objects;



    // Scene loading:
    std::string scene_Filename = "D:/RayTracing/RayTracing/rayTrace/rayTrace/data/scene.txt";
    ifstream textscene;
    textscene.open(scene_Filename);

    if (!textscene.is_open()) {
        exit(EXIT_FAILURE);
    }

    std::string word;
    std::string words[255];
    int wordcount = 0;

    textscene >> word;
    while (textscene.good()){
        if ((word[0] != '_')  && (word[1] != '_')) {
            words[wordcount].append(word);
         //   cout << words[wordcount] << " ";     // to check if every word is saved
        }
        textscene >> word;
        wordcount++;

    }
    int checker = 0;
    
    if (words[checker] == "Mesh") {
        checker++;
        if (words[checker] == "{") {
            checker++;
        }
        if (words[checker] == "yes") {
            meshIncluded = 1;
        }
        checker++;
        objFilename = words[checker];
        checker++;
        if (words[checker] == "}") {
            //cout << "Mesh File route: " << objFilename << endl;
            checker++;
        }
    }

    if (words[checker] == "GlobalSettings") {     // we read GlobalSettings
        checker++;

        if (words[checker] == "{") {
            checker++;
        }
        else { cout << "Error: Scene.txt GlobalSettings" << endl; }

        if (words[checker] == "aspect_ratio") {
            checker++;
            double a = std::stod(words[checker]);
            checker++;
            checker++;
            double b = std::stod(words[checker]);
            aspect_ratio = a/b; // stod() converts string to a double
            checker++;
        } // else aspect ratio not setted; maybe break;

        if (words[checker] == "image_width") {
            checker++;
            image_width = std::stoi(words[checker]); // stoi() makes int from a string
            checker++;
        } // else image_width is not setted; maybe break;

        if (words[checker] == "samples_per_pixel") {
            checker++;
            samples_per_pixel = std::stoi(words[checker]);
            checker++;
        } // else samples_per_pixel is not defined corectly; maybe break;

        if (words[checker] == "max_depth") {
            checker++;
            max_depth = std::stoi(words[checker]);
            checker++;
        } // else max_depth is not deined.; mby break;

        if (words[checker] == "background") {
            checker++;
            checker++;
            double R = std::stod(words[checker]);
            checker++;
            double G = std::stod(words[checker]);
            checker++;
            double B = std::stod(words[checker]);
            checker++;
            checker++;
            background.e[0] = R;
            background.e[1] = G;
            background.e[2] = B;
        }
        if (words[checker] == "}") {
            checker++;
            //cout << "GlobalSettings read succesfully." << endl;
        } // else global settings didnt close succesfully; break?
    } // end of GlobalSettings

    // If there is no mesh load camera 
    if (!meshIncluded) {
        if (words[checker] == "Camera") {
            checker++;
            if (words[checker] == "{") {
                checker++;
            }
            else { cout << "Error: Scene.txt GlobalSettings" << endl; }

            if (words[checker] == "lookfrom") { // lookfrom ( 13 2 3 )
                checker++;
                if (words[checker] == "(") {
                    checker++;
                    //point3 lookfrom;    // declare at top
                    lookfrom.e[0] = std::stod(words[checker]);
                    checker++;
                    lookfrom.e[1] = std::stod(words[checker]);
                    checker++;
                    lookfrom.e[2] = std::stod(words[checker]);
                    checker++;
                    if (words[checker] == ")") {
                        //cout << "Camera point3 lookfrom created successfully" << endl;
                        checker++;
                    }
                }
            }
            if (words[checker] == "lookat") { // lookat ( 0 0 0 )
                checker++;
                if (words[checker] == "(") {
                    checker++;
                //     point3 lookat;    // declare at top
                    lookat.e[0] = std::stod(words[checker]);
                    checker++;
                    lookat.e[1] = std::stod(words[checker]);
                    checker++;
                    lookat.e[2] = std::stod(words[checker]);
                    checker++;
                    if (words[checker] == ")") {
                      //  cout << "Camera point3 lookat created successfully" << endl;
                        checker++;
                    }
                }
            }
            if (words[checker] == "vup") {
                checker++;
                if (words[checker] == "(") {
                    checker++;
                    //  vec3 vup;    // declare at top
                    vup.e[0] = std::stod(words[checker]);
                    checker++;
                    vup.e[1] = std::stod(words[checker]);
                    checker++;
                    vup.e[2] = std::stod(words[checker]);
                    checker++;
                    if (words[checker] == ")") {
                      //  cout << "Camera vec3 vup created successfully" << endl;
                        checker++;
                    }
                }
            }
            if (words[checker] == "vfov") {
                checker++;
                vfov = std::stod(words[checker]); // stod() converts string to a double
                checker++;
            } // else vfov ratio not setted; maybe break;
            if (words[checker] == "aperture") {
                checker++;
                aperture = std::stod(words[checker]); // stod() converts string to a double
                checker++;
            } // else aperture ratio not setted; maybe break;
            if (words[checker] == "dist_to_focus") {
                checker++;
                dist_to_focus = std::stod(words[checker]); // stod() converts string to a double
                checker++;
            } // else dist_to_focus ratio not setted; maybe break;
            if (words[checker] == "}") {
                //cout << "Camera settet corectly" << endl;
                checker++;
            }
        }


        //shared_ptr<material> material;
        if (words[checker] == "Objects") {
            checker++;
            // now we are at {
            checker++; 
            while (words[checker] != "}") {
                //std::string name = words[checker];
                checker++;
                if (words[checker] == "lambertian") {

                    checker++;
                    // skip the (
                    checker++;
                    double red1 = std::stod(words[checker]);
                    checker++;
                    double green1 = std::stod(words[checker]);
                    checker++;
                    double blue1 = std::stod(words[checker]);
                    checker++;
                    // we are at ) now
                    checker++;

                    // get lambertian color

                    auto matherial = make_shared<lambertian>(color(red1, green1, blue1));

                    if (words[checker] == "sphere") {
                        checker++;
                        checker++; // to skip (
                        double X1 = std::stod(words[checker]);
                        checker++;
                        double Y1 = std::stod(words[checker]);
                        checker++;
                        double Z1 = std::stod(words[checker]);
                        checker++;
                        checker++; // to skip )
                        double R1 = std::stod(words[checker]);
                        checker++;

                        objects.add(make_shared<sphere>(point3(X1, Y1, Z1), R1, matherial));



                    }
                    if (words[checker] == "yz_rect") {
                        checker++;
                        checker++;
                        double a = std::stod(words[checker]);
                        checker++;
                        double b = std::stod(words[checker]);
                        checker++;
                        double c = std::stod(words[checker]);
                        checker++;
                        double d = std::stod(words[checker]);
                        checker++;
                        double e = std::stod(words[checker]);
                        checker++;
                        checker++;

                        objects.add(make_shared<yz_rect>(a, b, c, d, e, matherial));

                    }
                    if (words[checker] == "xz_rect") {
                        checker++;
                        checker++;
                        double a = std::stod(words[checker]);
                        checker++;
                        double b = std::stod(words[checker]);
                        checker++;
                        double c = std::stod(words[checker]);
                        checker++;
                        double d = std::stod(words[checker]);
                        checker++;
                        double e = std::stod(words[checker]);
                        checker++;
                        checker++;

                        objects.add(make_shared<xz_rect>(a, b, c, d, e, matherial));

                    }
                    if (words[checker] == "xy_rect") {
                        checker++;
                        checker++;
                        double a = std::stod(words[checker]);
                        checker++;
                        double b = std::stod(words[checker]);
                        checker++;
                        double c = std::stod(words[checker]);
                        checker++;
                        double d = std::stod(words[checker]);
                        checker++;
                        double e = std::stod(words[checker]);
                        checker++;
                        checker++;

                        objects.add(make_shared<xy_rect>(a, b, c, d, e, matherial));

                    }
                    if (words[checker] == "box") {
                        checker++;
                        checker++;
                        double p1X = std::stod(words[checker]);
                        checker++;
                        double p1Y = std::stod(words[checker]);
                        checker++;
                        double p1Z = std::stod(words[checker]);
                        checker++;
                        checker++;
                        checker++;
                        double p2X = std::stod(words[checker]);
                        checker++;
                        double p2Y = std::stod(words[checker]);
                        checker++;
                        double p2Z = std::stod(words[checker]);
                        checker++;

                        objects.add(make_shared<box>(point3(p1X, p1Y, p1Z), point3(p2X, p2Y, p2Z), matherial));

                    }

                    //objects.add(make_shared<sphere>(point3(0, -10, 0), 10, make_shared<lambertian>(checker)));

                }
                else if (words[checker] == "metal") {
                    checker++;
                    // skip the (
                    checker++;
                    double red1 = std::stod(words[checker]);
                    checker++;
                    double green1 = std::stod(words[checker]);
                    checker++;
                    double blue1 = std::stod(words[checker]);
                    checker++;
                    // we are at ) now
                    checker++;
                    double fuzz = std::stod(words[checker]);
                    checker++;

                    // get metal color
                    auto matherial = make_shared<metal>(color(red1, green1, blue1), fuzz);

                    if (words[checker] == "sphere") {
                        checker++;
                        checker++; // to skip (
                        double X1 = std::stod(words[checker]);
                        checker++;
                        double Y1 = std::stod(words[checker]);
                        checker++;
                        double Z1 = std::stod(words[checker]);
                        checker++;
                        checker++; // to skip )
                        double R1 = std::stod(words[checker]);
                        checker++;

                        objects.add(make_shared<sphere>(point3(X1, Y1, Z1), R1, matherial));

                    }
                    if (words[checker] == "yz_rect") {
                        checker++;
                        checker++;
                        double a = std::stod(words[checker]);
                        checker++;
                        double b = std::stod(words[checker]);
                        checker++;
                        double c = std::stod(words[checker]);
                        checker++;
                        double d = std::stod(words[checker]);
                        checker++;
                        double e = std::stod(words[checker]);
                        checker++;
                        checker++;

                        objects.add(make_shared<yz_rect>(a, b, c, d, e, matherial));

                    }
                    if (words[checker] == "xz_rect") {
                        checker++;
                        checker++;
                        double a = std::stod(words[checker]);
                        checker++;
                        double b = std::stod(words[checker]);
                        checker++;
                        double c = std::stod(words[checker]);
                        checker++;
                        double d = std::stod(words[checker]);
                        checker++;
                        double e = std::stod(words[checker]);
                        checker++;
                        checker++;

                        objects.add(make_shared<xz_rect>(a, b, c, d, e, matherial));

                    }
                    if (words[checker] == "xy_rect") {
                        checker++;
                        checker++;
                        double a = std::stod(words[checker]);
                        checker++;
                        double b = std::stod(words[checker]);
                        checker++;
                        double c = std::stod(words[checker]);
                        checker++;
                        double d = std::stod(words[checker]);
                        checker++;
                        double e = std::stod(words[checker]);
                        checker++;
                        checker++;

                        objects.add(make_shared<xy_rect>(a, b, c, d, e, matherial));

                    }
                    if (words[checker] == "box") {
                        checker++;
                        checker++;
                        double p1X = std::stod(words[checker]);
                        checker++;
                        double p1Y = std::stod(words[checker]);
                        checker++;
                        double p1Z = std::stod(words[checker]);
                        checker++;
                        checker++;
                        checker++;
                        double p2X = std::stod(words[checker]);
                        checker++;
                        double p2Y = std::stod(words[checker]);
                        checker++;
                        double p2Z = std::stod(words[checker]);
                        checker++;

                        objects.add(make_shared<box>(point3(p1X, p1Y, p1Z), point3(p2X, p2Y, p2Z), matherial));

                    }
                }
                else if (words[checker] == "dielectric") {
                    checker++;
                    double refraction_index = std::stod(words[checker]);
                    checker++;

                    auto matherial = make_shared<dielectric>(refraction_index);

                    if (words[checker] == "sphere") {
                        checker++;
                        checker++; // to skip (
                        double X1 = std::stod(words[checker]);
                        checker++;
                        double Y1 = std::stod(words[checker]);
                        checker++;
                        double Z1 = std::stod(words[checker]);
                        checker++;
                        checker++; // to skip )
                        double R1 = std::stod(words[checker]);
                        checker++;

                        objects.add(make_shared<sphere>(point3(X1, Y1, Z1), R1, matherial));



                    }
                    if (words[checker] == "yz_rect") {
                        checker++;
                        checker++;
                        double a = std::stod(words[checker]);
                        checker++;
                        double b = std::stod(words[checker]);
                        checker++;
                        double c = std::stod(words[checker]);
                        checker++;
                        double d = std::stod(words[checker]);
                        checker++;
                        double e = std::stod(words[checker]);
                        checker++;
                        checker++;

                        objects.add(make_shared<yz_rect>(a, b, c, d, e, matherial));

                    }
                    if (words[checker] == "xz_rect") {
                        checker++;
                        checker++;
                        double a = std::stod(words[checker]);
                        checker++;
                        double b = std::stod(words[checker]);
                        checker++;
                        double c = std::stod(words[checker]);
                        checker++;
                        double d = std::stod(words[checker]);
                        checker++;
                        double e = std::stod(words[checker]);
                        checker++;
                        checker++;

                        objects.add(make_shared<xz_rect>(a, b, c, d, e, matherial));

                    }
                    if (words[checker] == "xy_rect") {
                        checker++;
                        checker++;
                        double a = std::stod(words[checker]);
                        checker++;
                        double b = std::stod(words[checker]);
                        checker++;
                        double c = std::stod(words[checker]);
                        checker++;
                        double d = std::stod(words[checker]);
                        checker++;
                        double e = std::stod(words[checker]);
                        checker++;
                        checker++;

                        objects.add(make_shared<xy_rect>(a, b, c, d, e, matherial));

                    }
                    if (words[checker] == "box") {
                        checker++;
                        checker++;
                        double p1X = std::stod(words[checker]);
                        checker++;
                        double p1Y = std::stod(words[checker]);
                        checker++;
                        double p1Z = std::stod(words[checker]);
                        checker++;
                        checker++;
                        checker++;
                        double p2X = std::stod(words[checker]);
                        checker++;
                        double p2Y = std::stod(words[checker]);
                        checker++;
                        double p2Z = std::stod(words[checker]);
                        checker++;

                        objects.add(make_shared<box>(point3(p1X, p1Y, p1Z), point3(p2X, p2Y, p2Z), matherial));

                    }
                    // get dielectric refraction
                }
                else if (words[checker] == "diffuse_light") {
                    checker++;
                    // skip the (
                    checker++;
                    double red1 = std::stod(words[checker]);
                    checker++;
                    double green1 = std::stod(words[checker]);
                    checker++;
                    double blue1 = std::stod(words[checker]);
                    checker++;
                    // we are at ) now
                    checker++;

                    // get lambertian color

                    auto matherial = make_shared<diffuse_light>(color(red1, green1, blue1));

                    if (words[checker] == "sphere") {
                        checker++;
                        checker++; // to skip (
                        double X1 = std::stod(words[checker]);
                        checker++;
                        double Y1 = std::stod(words[checker]);
                        checker++;
                        double Z1 = std::stod(words[checker]);
                        checker++;
                        checker++; // to skip )
                        double R1 = std::stod(words[checker]);
                        checker++;

                        objects.add(make_shared<sphere>(point3(X1, Y1, Z1), R1, matherial));
                    }
                    if (words[checker] == "yz_rect") {
                        checker++;
                        checker++;
                        double a = std::stod(words[checker]);
                        checker++;
                        double b = std::stod(words[checker]);
                        checker++;
                        double c = std::stod(words[checker]);
                        checker++;
                        double d = std::stod(words[checker]);
                        checker++;
                        double e = std::stod(words[checker]);
                        checker++;
                        checker++;

                        objects.add(make_shared<yz_rect>(a, b, c, d, e, matherial));

                    }
                    if (words[checker] == "xz_rect") {
                        checker++;
                        checker++;
                        double a = std::stod(words[checker]);
                        checker++;
                        double b = std::stod(words[checker]);
                        checker++;
                        double c = std::stod(words[checker]);
                        checker++;
                        double d = std::stod(words[checker]);
                        checker++;
                        double e = std::stod(words[checker]);
                        checker++;
                        checker++;

                        objects.add(make_shared<xz_rect>(a, b, c, d, e, matherial));

                    }
                    if (words[checker] == "xy_rect") {
                        checker++;
                        checker++;
                        double a = std::stod(words[checker]);
                        checker++;
                        double b = std::stod(words[checker]);
                        checker++;
                        double c = std::stod(words[checker]);
                        checker++;
                        double d = std::stod(words[checker]);
                        checker++;
                        double e = std::stod(words[checker]);
                        checker++;
                        checker++;

                        objects.add(make_shared<xy_rect>(a, b, c, d, e, matherial));

                    }
                    if (words[checker] == "box") {
                        checker++;
                        checker++;
                        double p1X = std::stod(words[checker]);
                        checker++;
                        double p1Y = std::stod(words[checker]);
                        checker++;
                        double p1Z = std::stod(words[checker]);
                        checker++;
                        checker++;
                        checker++;
                        double p2X = std::stod(words[checker]);
                        checker++;
                        double p2Y = std::stod(words[checker]);
                        checker++;
                        double p2Z = std::stod(words[checker]);
                        checker++;

                        objects.add(make_shared<box>(point3(p1X, p1Y, p1Z), point3(p2X, p2Y, p2Z), matherial));

                    }
                }

            } // if "}" then we exit
            checker++;
        }
    }
    




    int image_height = static_cast<int>(image_width / aspect_ratio); // this is skipped and automatically setted

    
	auto res = cudaSetDevice(0);



    
 
    

    // OptiX Works here if there is a valid mesh loaded from scene
    if (meshIncluded) {
        RTPcontexttype contextType;
        RTPbuffertype bufferType;
        int runOnCpu = 1;
        if (runOnCpu) {
            contextType = RTP_CONTEXT_TYPE_CPU;
            bufferType = RTP_BUFFER_TYPE_HOST;
        }
        else {
            contextType = RTP_CONTEXT_TYPE_CUDA;
            bufferType = RTP_BUFFER_TYPE_CUDA_LINEAR;
        }

        RTPcontext context;
        if (runOnCpu) {
            CHK_PRIME(rtpContextCreate(contextType, &context));
        }
        else {
            CHK_PRIME(rtpContextCreate(contextType, &context));
            const unsigned deviceNumbers[] = { 0 };
            CHK_PRIME(rtpContextSetCudaDeviceNumbers(context, 1, deviceNumbers));
        }

        PrimeMesh mesh;
        loadMesh(objFilename, mesh);

        vec3 extents = make_vec3(mesh.getBBoxMax() - mesh.getBBoxMin());

        // Here we set camera automatically for Mesh scenes, becouse they are with different 
        // sizes and are hard to set from scene file.

        lookfrom = vec3(mesh.getBBoxMin().x, mesh.getBBoxMin().y, mesh.getBBoxMin().z);
        lookfrom = lookfrom - vec3(0.0, 0.0, -extents.length() * 2);
        lookat = vec3(mesh.getBBoxMin().x, mesh.getBBoxMin().y, mesh.getBBoxMin().z);
        vup = -unit_vector(cross(lookat - lookfrom, vec3(1, 0, 0)));
        auto dist_to_focus = (lookfrom - lookat).length();
        aperture = 0.0;
        vfov = 40.0;
        camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

        //
        // Create buffers for geometry data 
        //
        RTPbufferdesc indicesDesc;
        CHK_PRIME(rtpBufferDescCreate(
            context,
            RTP_BUFFER_FORMAT_INDICES_INT3,
            bufferType,
            mesh.getVertexIndices(),
            &indicesDesc)
        );
        CHK_PRIME(rtpBufferDescSetRange(indicesDesc, 0, mesh.num_triangles));


        RTPbufferdesc verticesDesc;
        CHK_PRIME(rtpBufferDescCreate(
            context,
            RTP_BUFFER_FORMAT_VERTEX_FLOAT3,
            bufferType,
            mesh.getVertexData(),
            &verticesDesc)
        );
        CHK_PRIME(rtpBufferDescSetRange(verticesDesc, 0, mesh.num_vertices));

        //
        // Create the Model object
        //
        RTPmodel model;
        CHK_PRIME(rtpModelCreate(context, &model));
        CHK_PRIME(rtpModelSetTriangles(model, indicesDesc, verticesDesc));
        CHK_PRIME(rtpModelUpdate(model, 0));

        //
        // Create buffer for ray input 
        //
        RTPbufferdesc raysDesc;
        Buffer<Ray> raysBuffer(size_t(image_width) * image_height * samples_per_pixel, bufferType, UNLOCKED); 	// unlocked here

        // vector ot vsi4ki rays
        int index = 0;
        for (int j = image_height - 1; j >= 0; --j) {
            for (int i = 0; i < image_width; ++i) {
                color pixel_color(0, 0, 0);
                for (int s = 0; s < samples_per_pixel; ++s) {
                    auto u = (i + random_double()) / (image_width - 1);
                    auto v = (j + random_double()) / (image_height - 1);
                    ray r = cam.get_ray(u, v);
                    Ray Ray;

                    Ray.dir = make_float3(unit_vector(r.dir));
                    Ray.origin = make_float3(r.orig);
                    Ray.tmin = 0;
                    Ray.tmax = 1e+10;
                    raysBuffer.ptr()[index] = Ray;
                    index++;
                    //pixel_color += rt_color(h, Ray, background, mesh);
                }
            }
        }



        CHK_PRIME(rtpBufferDescCreate(
            context,
            Ray::format, /*RTP_BUFFER_FORMAT_RAY_ORIGIN_TMIN_DIRECTION_TMAX*/
            raysBuffer.type(),
            raysBuffer.ptr(),
            &raysDesc)
        );
        CHK_PRIME(rtpBufferDescSetRange(raysDesc, 0, raysBuffer.count()));


        //
        // Create buffer for returned hit descriptions
        //
        RTPbufferdesc hitsDesc;
        Buffer<Hit> hitsBuffer(raysBuffer.count(), bufferType, UNLOCKED);
        CHK_PRIME(rtpBufferDescCreate(
            context,
            Hit::format, /*RTP_BUFFER_FORMAT_HIT_T_TRIID_U_V*/
            hitsBuffer.type(),
            hitsBuffer.ptr(),
            &hitsDesc)
        );
        CHK_PRIME(rtpBufferDescSetRange(hitsDesc, 0, hitsBuffer.count()));


        //
        // Execute query
        //
        RTPquery query;
        CHK_PRIME(rtpQueryCreate(model, RTP_QUERY_TYPE_CLOSEST, &query));
        CHK_PRIME(rtpQuerySetRays(query, raysDesc));
        CHK_PRIME(rtpQuerySetHits(query, hitsDesc));
        CHK_PRIME(rtpQueryExecute(query, 0 /* hints */));

        //background = color(0.70, 0.80, 1.00);
        //std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";


        for (int i = 0; i < hitsBuffer.count(); i++) {
            Hit h = hitsBuffer.ptr()[i];
            Ray r = raysBuffer.ptr()[i];
            // color pixel_color(0.0, 0.0, 0.0);

           //  pixel_color = rt_color(h, r, background, mesh);
           //  write_rt_color(std::cout, pixel_color);

        }

        // Render
        //std::cout << "Writing results to file image123.ppm" << std::endl;

        // Outpitting the result directly into ppm file
        std::ofstream file("image123.ppm", std::ios::out);
        file << "P3\n" << image_width << ' ' << image_height << "\n255\n";
        int asd = 0;
        for (int j = image_height - 1; j >= 0; --j) {
            //    std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
            for (int i = 0; i < image_width; ++i) {
                color pixel_color(0, 0, 0);
                for (int s = 0; s < samples_per_pixel; ++s) {
                    const Hit& h = hitsBuffer.ptr()[asd++];
                    if (h.t > 0) {
                        int3* triangles = mesh.getVertexIndices();
                        int iv0 = triangles[h.triId].x;
                        int iv1 = triangles[h.triId].y;
                        int iv2 = triangles[h.triId].z;

                        float3* verts = mesh.getVertexData();
                        float3 A = verts[iv0];
                        float3 B = verts[iv1];
                        float3 C = verts[iv2];

                        vec3 AC = make_vec3(C - A);
                        vec3 AB = make_vec3(B - A);

                        vec3 n = unit_vector(cross(AC, AB));
                        //float3 N = normalize(::cross(AB,AC));
                        pixel_color[0] = pixel_color[0] + abs(n[0]);
                        pixel_color[1] = pixel_color[1] + abs(n[1]);
                        pixel_color[2] = pixel_color[2] + abs(n[2]);
                    }
                }
                pixel_color /= samples_per_pixel;
                write_color(file, pixel_color);
            }
        }
        file.close();

        freeMesh(mesh);

        //
        // cleanup
        //
        CHK_PRIME(rtpContextDestroy(context));
        // at bottom! 

    }



    // Base render if there isn't mesh loaded from scene
    if (!meshIncluded) {

        camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);


        //std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";
        std::ofstream file("image234.ppm", std::ios::out);
        file << "P3\n" << image_width << ' ' << image_height << "\n255\n";

        for (int j = image_height - 1; j >= 0; --j) {
            //    std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
            for (int i = 0; i < image_width; ++i) {
                color pixel_color(0, 0, 0);
                for (int s = 0; s < samples_per_pixel; ++s) {
                    auto u = (i + random_double()) / (image_width - 1);
                    auto v = (j + random_double()) / (image_height - 1);
                    ray r = cam.get_ray(u, v);
                    pixel_color += ray_color(r, background, objects, max_depth);
                }
                write_color(file, pixel_color, samples_per_pixel);
            }
        }
        file.close();



    
    }
    
}