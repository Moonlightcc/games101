#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.
        auto incrument = (end-start)/(num_nodes*1.0f-1.0f);
        
        for(int i = 0; i < num_nodes; ++i)
            masses.push_back(new Mass(start+i*incrument,node_mass,false));
        for(int i = 0; i < num_nodes-1; ++i)
            springs.push_back(new Spring(masses[i],masses[i+1],k));
        for (auto &i : pinned_nodes) {
            masses[i]->pinned = true;
        }
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            auto dis = s->m1->position - s->m2->position;
            auto mag = dis.norm();

            s->m1->forces += -dis/mag*(mag - s->rest_length) * s->k;
            s->m2->forces += dis/mag*(mag - s->rest_length) * s->k;
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position
                m->forces += m->mass * gravity;
                // TODO (Part 2): Add global damping
                float k_d_global = 0.01;
                m->forces += - k_d_global * m->velocity;

                auto a = m->forces / m->mass;
                m->velocity += a * delta_t;
                m->position += m->velocity * delta_t;
            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            auto dis = s->m1->position - s->m2->position;
            auto mag = dis.norm();

            s->m1->forces += -dis/mag*(mag - s->rest_length) * s->k;
            s->m2->forces += dis/mag*(mag - s->rest_length) * s->k;
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                m->forces += gravity * m->mass;
                Vector2D a = m->forces / m->mass;

                // TODO (Part 3.1): Set the new position of the rope mass
                Vector2D lastposition = m->position;
                // TODO (Part 4): Add global Verlet damping
                float dampfactor = 0.00005;
                m->position = m->position +  (1 - dampfactor) * (m->position - m->last_position) + a * delta_t *delta_t;
                m->last_position = lastposition;
            }
            m->forces =  Vector2D(0,0);
        }
    }
}
