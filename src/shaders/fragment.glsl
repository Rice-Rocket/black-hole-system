#version 410

#define MAX_STEPS 40
#define MAX_DIST 256.0
#define MIN_DIST 0.01
#define EPSILON 0.01
#define PI 3.1415927

#define SPEED 0.5
#define STEPS 12.0
#define SIZE 0.3

#define MIN_STEP_SIZE 0.375
#define STEP_CURVATURE_FACTOR 4.0

in vec2 fragCoord;
out vec4 fragColor;

uniform float time;
uniform vec2 resolution;
uniform vec4 mouse;

uniform mat3 camera;
uniform vec3 camera_origin;
uniform float camera_focal_length;

uniform sampler2D starmap;
uniform sampler2D rgba_noise;
uniform sampler2D organic_tex;

mat4 rotate_mat(vec3 theta) {
    float yaw = theta.x;
    float pitch = theta.y;
    float roll = theta.z;

    float cos_yaw = cos(yaw);
    float sin_yaw = sin(yaw);
    float cos_pitch = cos(pitch);
    float sin_pitch = sin(pitch);
    float cos_roll = cos(roll);
    float sin_roll = sin(roll);

    return mat4(
        vec4(cos_yaw, -sin_yaw, 0.0, 0.0),
        vec4(sin_yaw, cos_yaw, 0.0, 0.0),
        vec4(0.0, 0.0, 1.0, 0.0),
        vec4(0.0, 0.0, 0.0, 1.0)
    ) * mat4(
        vec4(cos_pitch, 0.0, sin_pitch, 0.0),
        vec4(0.0, 1.0, 0.0, 0.0),
        vec4(-sin_pitch, 0.0, cos_pitch, 0.0),
        vec4(0.0, 0.0, 0.0, 1.0)
    ) * mat4(
        vec4(1.0, 0.0, 0.0, 0.0),
        vec4(0.0, cos_roll, -sin_roll, 0.0),
        vec4(0.0, sin_roll, cos_roll, 0.0),
        vec4(0.0, 0.0, 0.0, 1.0)
    );
}

vec3 transform(vec3 p, vec3 theta, vec3 offset) {
    mat4 rot = rotate_mat(theta);
    return (rot * vec4(p - offset, 1.0)).xyz;
}

float hash(float x) {
    return fract(sin(x) * 152754.742);
}

float hash(vec2 x) {
    return hash(x.x + hash(x.y));
}

float pcurve( float x, float a, float b ) {
    float k = pow(a+b,a+b) / (pow(a,a)*pow(b,b));
    return k * pow( x, a ) * pow( 1.0-x, b );
}

float atan2(float y, float x)
{
    if (x > 0.0)
    {
        return atan(y / x);
    }
    else if (x == 0.0)
    {
        if (y > 0.0)
        {
            return PI / 2.0;
        }
        else if (y < 0.0)
        {
            return -(PI / 2.0);
        }
        else
        {
            return 0.0;
        }
    }
    else //(x < 0.0)
    {
        if (y >= 0.0)
        {
            return atan(y / x) + PI;
        }
        else
        {
            return atan(y / x) - PI;
        }
    }
}

float length_squared(vec3 x) {
    return dot(x, x);
}

float valnoise(vec2 p, float f) {
    float bl = hash(floor(p * f + vec2(0.0, 0.0)));
    float br = hash(floor(p * f + vec2(1.0, 0.0)));
    float tl = hash(floor(p * f + vec2(0.0, 1.0)));
    float tr = hash(floor(p * f + vec2(1.0, 1.0)));

    vec2 fr = fract(p * f);
    fr = (3.0 - 2.0 * fr) * fr * fr;
    float b = mix(bl, br, fr.x);
    float t = mix(tl, tr, fr.x);
    return mix(b, t, fr.y);
}

float valnoise(vec3 p) {
    float bl = hash(floor(p.xy * p.z + vec2(0.0, 0.0)));
    float br = hash(floor(p.xy * p.z + vec2(1.0, 0.0)));
    float tl = hash(floor(p.xy * p.z + vec2(0.0, 1.0)));
    float tr = hash(floor(p.xy * p.z + vec2(1.0, 1.0)));

    vec2 fr = fract(p.xy * p.z);
    fr = (3.0 - 2.0 * fr) * fr * fr;
    float b = mix(bl, br, fr.x);
    float t = mix(tl, tr, fr.x);
    return mix(b, t, fr.y);
}

float noise(in vec3 x) {
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f * f * (3.0 - 2.0 * f);
    vec2 uv = (p.xy + vec2(37.0,17.0) * p.z) + f.xy;
    vec2 rg = textureLod(rgba_noise, fract((uv+ 0.5) / 256.0), 0.0).yx;
    return -1.0 + 2.0 * mix(rg.x, rg.y, f.z);
}

vec4 background(vec3 ray) {
    vec2 uv = ray.xy;
    if (abs(ray.x) > 0.5) uv.x = ray.z;
    else if (abs(ray.y) > 0.5) uv.y = ray.z;

    float brightness = valnoise(uv * 3.0, 100.0);
    float color = valnoise(uv * 2.0, 20.0);

    brightness = pow(brightness, 256.0);
    brightness = brightness * 100.0;
    brightness = clamp(brightness, 0.0, 1.0);

    vec3 stars = brightness * mix(vec3(1.0, 0.6, 0.2), vec3(0.2, 0.6, 1.0), color);
    vec4 nebulae = texture(starmap, fract(uv * 1.5));
    nebulae.xyz += nebulae.xxx + nebulae.yyy + nebulae.zzz;
    nebulae.xyz *= 0.25;

    nebulae *= nebulae;
    nebulae *= nebulae;
    nebulae *= nebulae;
    nebulae *= nebulae;

    nebulae.xyz += stars;
    return nebulae;
}

float sd_torus(vec3 p, vec2 t) {
    vec2 q = vec2(length(p.xz) - t.x, p.y);
    return length(q)-t.y;
}

float sd_sphere(vec3 p, float r) {
  return length(p) - r;
}

vec2 sphere_intersect(vec3 ro, vec3 rd, vec3 center, float radius) {
    vec3 oc = ro - center;
    float a = dot(rd, rd);
    float half_b = dot(oc, rd);
    float c = dot(oc, oc) - radius * radius;

    float discriminant = half_b * half_b - a * c;
    if (discriminant < 0.0) return vec2(0.0);
    float sqrtd = sqrt(discriminant);

    float root = (-half_b - sqrtd) / a;
    // if (root < t_min || t_max < root) {
    //     root = (-half_b + sqrtd) / a;
    //     if (root < t_min || t_max < root) return vec2(0.0);
    // }
    return vec2(1.0, root);
}


void haze(inout vec3 color, float alpha, vec3 pos) {
    vec2 t = vec2(1.0, 0.01);

    float torus_dist = length(sd_torus(pos + vec3(0.0, -0.05, 0.0), t));

    float bloom_disc = 1.0 / (pow(torus_dist, 2.0) + 0.001);
    vec3 col = vec3(1.0, 1.0, 1.0);
    bloom_disc *= length(pos) < 0.5 ? 0.0 : 1.0;

    color += 0.01 * col * bloom_disc * (2.9 / float(MAX_STEPS)) * (1.0 - alpha * 1.0);
}

void accretion_disc(inout vec3 color, inout float alpha, vec3 pos) {
    float disc_radius = 3.2;
    float disc_width = 5.3;
    float disc_inner = disc_radius - disc_width * 0.5;
    float disc_outer = disc_radius + disc_width * 0.5;
    
    vec3 origin = vec3(0.0, 0.0, 0.0);
    vec3 disc_normal = normalize(vec3(0.0, 1.0, 0.0));
    float disc_thickness = 0.1;

    float center_dist = distance(pos, origin);
    float disc_dist = dot(disc_normal, pos - origin);
    
    float rad_grad = 1.0 - clamp((center_dist - disc_inner) / disc_width * 0.5, 0.0, 1.0);

    float coverage = pcurve(rad_grad, 4.0, 0.9);

    disc_thickness *= rad_grad;
    coverage *= clamp(1.0 - abs(disc_dist) / disc_thickness, 0.0, 1.0);

    vec3 dust_lit = vec3(1.0, 1.0, 1.0);
    vec3 dust_dark = vec3(0.0, 0.0, 0.0);

    float dust_glow = 1.0 / (pow(1.0 - rad_grad, 2.0) * 290.0 + 0.002);
    vec3 dust_col = dust_lit * dust_glow * 8.2;

    coverage = clamp(coverage * 0.7, 0.0, 1.0);


    float fade = pow((abs(center_dist - disc_inner) + 0.4), 4.0) * 0.04;
    float bloom_factor = 1.0 / (pow(disc_dist, 2.0) * 40.0 + fade + 0.00002);
    vec3 b = dust_lit * pow(bloom_factor, 1.5);
    
    b *= mix(vec3(1.7, 1.1, 1.0), vec3(0.5, 0.6, 1.0), vec3(pow(rad_grad, 2.0)));
    b *= mix(vec3(1.7, 0.5, 0.1), vec3(1.0), vec3(pow(rad_grad, 0.5)));

    dust_col = mix(dust_col, b * 150.0, clamp(1.0 - coverage * 1.0, 0.0, 1.0));
    coverage = clamp(coverage + bloom_factor * bloom_factor * 0.1, 0.0, 1.0);
    
    if (coverage < 0.01) return;   
    
    vec3 rad_coords;
    rad_coords.x = center_dist * 1.5 + 0.55;
    rad_coords.y = atan2(-pos.x, -pos.z) * 1.5;
    rad_coords.z = disc_dist * 1.5;

    rad_coords *= 0.95;
    float speed = 1.0;
    
    float n1 = 1.0;
    vec3 rc = rad_coords + 0.0;             rc.y += time * speed;
    n1 *= noise(rc * 3.0) * 0.5 + 0.5;      rc.y -= time * speed;
    n1 *= noise(rc * 6.0) * 0.5 + 0.5;      rc.y += time * speed;
    n1 *= noise(rc * 12.0) * 0.5 + 0.5;     rc.y -= time * speed;
    n1 *= noise(rc * 24.0) * 0.5 + 0.5;     rc.y += time * speed;

    float n2 = 2.0;
    rc = rad_coords + 30.0;
    n2 *= noise(rc * 3.0) * 0.5 + 0.5;      rc.y += time * speed;
    n2 *= noise(rc * 6.0) * 0.5 + 0.5;      rc.y -= time * speed;
    n2 *= noise(rc * 12.0) * 0.5 + 0.5;     rc.y += time * speed;
    n2 *= noise(rc * 24.0) * 0.5 + 0.5;     rc.y -= time * speed;
    n2 *= noise(rc * 48.0) * 0.5 + 0.5;     rc.y += time * speed;
    n2 *= noise(rc * 92.0) * 0.5 + 0.5;     rc.y -= time * speed;

    dust_col *= n1 * 0.998 + 0.002;
    coverage *= n2;
    
    rad_coords.y += time * speed * 0.5;
    dust_col *= pow(texture(organic_tex, fract(rad_coords.yx * vec2(0.15, 0.27))).rgb, vec3(2.0)) * 4.0;

    coverage = clamp(coverage * 1200.0 / float(MAX_STEPS), 0.0, 1.0);
    dust_col = max(vec3(0.0), dust_col);
    coverage *= pcurve(rad_grad, 4.0, 0.9);

    color = (1.0 - alpha) * dust_col * coverage + color;
    alpha = (1.0 - alpha) * coverage + alpha;
}


vec4 scene(vec3 p) {
    vec4 res = vec4(1.0, 0.5, 0.25, sd_sphere(p - vec3(7.0, 1.0, 3.0), 0.2));
    return res;
}


vec4 raymarch_scene(vec3 ro, vec3 rd) {
    float dist = 0;
    vec4 res = vec4(-1.0);

    for (int i = 0; i < MAX_STEPS && dist <= MAX_DIST; i++) {
        vec3 p = ro + rd * dist;
        vec4 ds = scene(p);
        if (abs(ds.w) < (MIN_DIST * dist)) {
            res.w = dist;
            res.rgb = ds.rgb;
            break;
        }
        dist += ds.w;
    }
    if (res.w < -0.5) res.rgb = background(rd).rgb;

    return res;
}

vec4 render(vec3 rd, vec2 uv) {
    const float far = 15.0;
    vec3 bh_origin = vec3(0.0, 0.0, 0.0);
    float bh_bound_radius = 10.0;

    vec3 color = vec3(0.0, 0.0, 0.0);
    float alpha = 0.0;
    float dither = clamp(fract(sin(dot(uv, vec2(12.9898, 78.223))) * 43758.5453), 0.0, 1.0) * 1.0;

    float stepsize = MIN_STEP_SIZE;
    float curvature;
    float dist = 0.0;
    float resdist = -1.0;
    vec3 raypos = camera_origin + rd * dither * MIN_STEP_SIZE;
    // raypos += rd * dither * sd_sphere(raypos - bh_origin, bh_bound_radius);
    for (int i = 0; i < MAX_STEPS; i++) {
        // gravitational lensing
        float singularity_dist = distance(raypos, bh_origin);
        float warp_factor = 1.0 / (pow(singularity_dist, 2.0) + 0.000001);
        vec3 singularity_vec = normalize(bh_origin - raypos);
        float warp_amount = 5.0;
        curvature = length(rd - singularity_vec);
        rd = normalize(rd + singularity_vec * warp_factor * warp_amount / float(MAX_STEPS));
        if (singularity_dist < bh_bound_radius) stepsize = MIN_STEP_SIZE;
        else stepsize = max(singularity_dist / STEP_CURVATURE_FACTOR, MIN_STEP_SIZE);
        // if (length(raypos - bh_origin) > bh_bound_radius) break;
        
        // accretion disc
        accretion_disc(color, alpha, raypos);
        haze(color, alpha, raypos);

        // scene
        vec4 scene_dist = scene(raypos);
        if (abs(scene_dist.w) < (MIN_DIST * dist)) {
            color += scene_dist.rgb;
            resdist = dist;
            break;
        }
        dist += scene_dist.w;
        raypos += rd * min(scene_dist.w, stepsize);
    }
    // color *= 1.0 / float(MAX_STEPS) / 1.0;
    if (resdist < -0.5 && (alpha < 0.1 && length(raypos - bh_origin) > 0.5)) color += background(rd).rgb;

    color = pow(color, vec3(1.5));
    color = color / (1.0 + color);
    color = pow(color, vec3(1.0 / 1.5));
    color = mix(color, color * color * (3.0 - 2.0 * color), vec3(1.0));
    color = pow(color, vec3(1.3, 1.2, 1.0));
    color = clamp(color * 1.01, 0.0, 1.0);
    color = pow(color, vec3(0.7 / 2.2));

    // color += (alpha > 0.1 || length(raypos - bh_origin) < 0.5) ? vec3(0.0) : raymarch_scene(raypos, rd).rgb;
    return vec4(clamp(color, 0.0, 1.0), 1.0);
}

void main() {
    vec2 uv = ((fragCoord * resolution.xy) - 0.5 * resolution.xy) / resolution.y;
    vec3 ray_direction = camera * normalize(vec3(uv, camera_focal_length));
    vec4 color = render(ray_direction, uv);
    fragColor = color;
}