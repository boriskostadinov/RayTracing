//#include "color.h"
//#include "ray.h"
//#include "vec3.h"


#include<fstream>

#include "rtweekend.h"
#include "camera.h"
#include "color.h"
#include "hittable_list.h"
#include "material.h"
#include "moving_sphere.h"
#include "sphere.h"
#include "aarect.h"

#include <iostream>

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
	auto res = cudaSetDevice(0);
    // Image 
/*
const auto aspect_ratio = 16.0 / 9.0;
const int image_width = 400;
const int image_height = static_cast<int>(image_width / aspect_ratio);
const int samples_per_pixel = 100;
const int max_depth = 50;
*/
    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 800;
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 10;
    const int max_depth = 50;



    // World
    /*hittable_list world;

    auto material_ground = make_shared<lambertian>(color(0.8, 0.8, 0.0));
    auto material_center = make_shared<lambertian>(color(0.1, 0.2, 0.5));
    auto material_left = make_shared<dielectric>(1.5);
    auto material_right = make_shared<metal>(color(0.8, 0.6, 0.2), 0.0);


    world.add(make_shared<sphere>(point3(0.0, -100.5, -1.0), 100.0, material_ground));
    world.add(make_shared<sphere>(point3(0.0, 0.0, -1.0), 0.5, material_center));
    world.add(make_shared<sphere>(point3(-1.0, 0.0, -1.0), 0.5, material_left));
    world.add(make_shared<sphere>(point3(-1.0, 0.0, -1.0), -0.4, material_left));
    world.add(make_shared<sphere>(point3(1.0, 0.0, -1.0), 0.5, material_right));

    //world.add(make_shared<sphere>(point3(0, 0, -1), 0.5));
    //world.add(make_shared<sphere>(point3(0, -100.5, -1), 100));
    */
    /*
    auto R = cos(pi / 4);
    hittable_list world;

    auto material_left = make_shared<lambertian>(color(0, 0, 1));
    auto material_right = make_shared<lambertian>(color(1, 0, 0));

    world.add(make_shared<sphere>(point3(-R, 0, -1), R, material_left));
    world.add(make_shared<sphere>(point3(R, 0, -1), R, material_right));
    */
    /*
    hittable_list world;

    auto material_ground = make_shared<lambertian>(color(0.8, 0.8, 0.0));
    auto material_center = make_shared<lambertian>(color(0.1, 0.2, 0.5));
    auto material_left = make_shared<dielectric>(1.5);
    auto material_right = make_shared<metal>(color(0.8, 0.6, 0.2), 0.0);

    world.add(make_shared<sphere>(point3(0.0, -100.5, -1.0), 100.0, material_ground));
    world.add(make_shared<sphere>(point3(0.0, 0.0, -1.0), 0.5, material_center));
    world.add(make_shared<sphere>(point3(-1.0, 0.0, -1.0), 0.5, material_left));
    world.add(make_shared<sphere>(point3(-1.0, 0.0, -1.0), -0.45, material_left));
    world.add(make_shared<sphere>(point3(1.0, 0.0, -1.0), 0.5, material_right));

    //camera cam(point3(-2, 2, 1), point3(0, 0, -1), vec3(0, 1, 0), 90, aspect_ratio);
    //camera cam(point3(-2, 2, 1), point3(0, 0, -1), vec3(0, 1, 0), 20, aspect_ratio);

    point3 lookfrom(3, 3, 2);
    point3 lookat(0, 0, -1);
    vec3 vup(0, 1, 0);
    auto dist_to_focus = (lookfrom - lookat).length();
    auto aperture = 2.0;

    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus);
    */
    hittable_list world;

    point3 lookfrom;
    point3 lookat;
    auto vfov = 40.0;
    auto aperture = 0.0;
    color background(0, 0, 0);


    switch (0) {
    case 1:
        world = random_scene();
        background = color(0.70, 0.80, 1.00);
        lookfrom = point3(13, 2, 3);
        lookat = point3(0, 0, 0);
        vfov = 20.0;
        aperture = 0.1;
        break;

    case 2:
        world = two_spheres();
        background = color(0.70, 0.80, 1.00);
        lookfrom = point3(13, 2, 3);
        lookat = point3(0, 0, 0);
        vfov = 20.0;
        break;

    case 4:
        world = earth();
        background = color(0.70, 0.80, 1.00);
        lookfrom = point3(13, 2, 3);
        lookat = point3(0, 0, 0);
        vfov = 20.0;
        break;

    case 5:
        world = simple_light();
        //samples_per_pixel = 400;
        background = color(0, 0, 0);
        lookfrom = point3(26, 3, 6);
        lookat = point3(0, 2, 0);
        vfov = 20.0;
        break;
    default:
        break;

    }

    // Camera
    /*auto viewport_height = 2.0;
    auto viewport_width = aspect_ratio * viewport_height;
    auto focal_length = 1.0;

    auto origin = point3(0, 0, 0);
    auto horizontal = vec3(viewport_width, 0, 0);
    auto vertical = vec3(0, viewport_height, 0);
    auto lower_left_corner = origin - horizontal / 2 - vertical / 2 - vec3(0, 0, focal_length);
    */
    /*point3 lookfrom(13, 2, 3);
    point3 lookat(0, 0, 0);
    vec3 vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.1;
    int image_height = static_cast<int>(image_width / aspect_ratio);

    camera cam(lookfrom, lookat, vup, 20, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);
*/
    vec3 vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    //int image_height = static_cast<int>(image_width / aspect_ratio);


	RTPcontexttype contextType;
	RTPbuffertype bufferType;
	int runOnCpu = 1;
	if (runOnCpu) {
		contextType = RTP_CONTEXT_TYPE_CPU;
		bufferType = RTP_BUFFER_TYPE_HOST;
	} else {
		contextType = RTP_CONTEXT_TYPE_CUDA;
		bufferType = RTP_BUFFER_TYPE_CUDA_LINEAR;
	}

	RTPcontext context;
	if (runOnCpu) {
		CHK_PRIME(rtpContextCreate(contextType, &context));
	} else {
		CHK_PRIME(rtpContextCreate(contextType, &context));
		const unsigned deviceNumbers[] = {0};
		CHK_PRIME(rtpContextSetCudaDeviceNumbers(context, 1, deviceNumbers));
	}

    std::string objFilename = "D:/RayTracing/RayTracing/rayTrace/rayTrace/data/cow.obj";
    //std::string objFilename = "/data/cow.obj";

    PrimeMesh mesh;
    loadMesh(objFilename, mesh);

	vec3 extents = make_vec3(mesh.getBBoxMax() - mesh.getBBoxMin());

	lookfrom = vec3(mesh.getBBoxMin().x, mesh.getBBoxMin().y, mesh.getBBoxMin().z);
	lookfrom = lookfrom - vec3(0.0,0.0, -extents.length()*2);
    lookat = vec3(mesh.getBBoxMin().x,mesh.getBBoxMin().y,mesh.getBBoxMin().z);
	vup = -unit_vector(cross(lookat-lookfrom, vec3(1,0,0)));
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
    Buffer<Ray> raysBuffer(size_t(image_width)*image_height*samples_per_pixel, bufferType, UNLOCKED); 	// unlocked here
    
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

                Ray.dir    = make_float3(unit_vector(r.dir));
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
    std::cout << "Writing results to file image123.ppm" << std::endl;

    // Outpitting the result directly into ppm file
    std::ofstream file("image123.ppm", std::ios::out);
    file << "P3\n" << image_width << ' ' << image_height << "\n255\n";
	int asd = 0;
	for (int j = image_height - 1; j >= 0; --j) {
        //    std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
				const Hit &h = hitsBuffer.ptr()[asd++];
				if (h.t > 0) {
					int3* triangles = mesh.getVertexIndices();
					int iv0 = triangles[h.triId].x;
					int iv1 = triangles[h.triId].y;
					int iv2 = triangles[h.triId].z;

					float3* verts = mesh.getVertexData();
					float3 A = verts[iv0];
					float3 B = verts[iv1];
					float3 C = verts[iv2];

					vec3 AC = make_vec3(C-A);
					vec3 AB = make_vec3(B-A);
					
					vec3 n = unit_vector(cross(AC,AB));
					//float3 N = normalize(::cross(AB,AC));
					pixel_color[0] = pixel_color[0]+abs(n[0]);
					pixel_color[1] = pixel_color[1]+abs(n[1]);
					pixel_color[2] = pixel_color[2]+abs(n[2]);
				}
			}
			pixel_color /= samples_per_pixel;
			write_color(file, pixel_color);
		}
	}
    file.close();

    // Not included TODO: my render fx();
   
    //
    // Shade the hit results to create image
    //
   
    //std::vector<float3> image(image_width * image_height);
    //shadeHits(image, hitsBuffer, mesh);
    //writePpm("output.ppm", &image[0].x, image_width, image_height);

// HERE:





    //
    // re-execute query with different rays
    //
    //translateRays(raysBuffer, extents * make_float3(0.2f, 0, 0));
    //CHK_PRIME(rtpQueryExecute(query, 0 /* hints */));
    //shadeHits(image, hitsBuffer, mesh);
    freeMesh(mesh);
    //writePpm("outputTranslated.ppm", &image[0].x, image_width, image_height);




    //
    // cleanup
    //
    CHK_PRIME(rtpContextDestroy(context));
    // at bottom! 
    





    // Render
    /*
    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";

    // Outpitting the result directly into ppm file
    //std::ofstream file("image123.ppm", std::ios::out);
    //file << "P3\n" << image_width << ' ' << image_height << "\n255\n";
    //write_color(std::cout, pixel_color);  must be into file << pixel_color;
    //file.close();

    for (int j = image_height - 1; j >= 0; --j) {
        //    std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width - 1);
                auto v = (j + random_double()) / (image_height - 1);
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, background, world, max_depth);
            }
            write_color(std::cout, pixel_color, samples_per_pixel);
        }
    }

    // std::cerr << "\nDone.\n";
    */
}