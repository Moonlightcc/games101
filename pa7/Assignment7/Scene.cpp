//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, const Intersection &interp, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    if(interp.m->hasEmission())
        return interp.m->getEmission();
    
    auto shade_point = interp.coords;
    auto wo = -ray.direction;
    auto mat = interp.m;
    auto normal = normalize(interp.normal);

    Intersection light_point;
    float pdf;
    sampleLight(light_point,pdf);
    auto sample_point = light_point.coords;
    auto shadep_to_lightp = sample_point - shade_point;
    auto ws = normalize(shadep_to_lightp);
    auto light_normal = normalize(light_point.normal);

    Vector3f L_dir;
    auto inter_dir = intersect(Ray(shade_point, ws));
    float dist_dir = (inter_dir.coords - light_point.coords).norm();
    if ( dist_dir < 0.0005f )
    {
        L_dir = light_point.emit * mat->eval(ws, wo, normal) * std::max(0.0f,dotProduct(ws, normal)) 
        * std::max(0.0f,dotProduct(-ws, light_normal)) / dotProduct(shadep_to_lightp, shadep_to_lightp) / pdf;
    }

    Vector3f L_indir{};
    if(get_random_float()<RussianRoulette)
    {
        auto wi = mat->sample(wo, normal).normalized();
        auto pdf = mat->pdf(wo, wi, normal);
        if(pdf > 0.0005f)
        {
            auto new_ray = Ray(shade_point, wi);
            auto nextp = intersect(new_ray);
            if (nextp.happened && !nextp.m->hasEmission())
            {
                auto brdf = mat->eval(wi, wo, normal);
                L_indir = castRay(new_ray,nextp,depth) * brdf * std::max(.0f, dotProduct(wi, normal)) / pdf / RussianRoulette;
            }
        }
    }
    return L_dir + L_indir;
}