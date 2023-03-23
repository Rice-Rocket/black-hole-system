#version 410

#define MAX_STEPS 64
#define MAX_DIST 10.0
#define MIN_DIST 0.01
#define EPSILON 0.01
#define PI 3.1415926

#define SPEED 0.5
#define STEPS 12.0
#define SIZE 0.3

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
// uniform sampler2D organic_tex;

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

vec4 raymarch_disk(vec3 ray, vec3 pos) {
    float disc_radius = 0.6;
    float disc_width = 1.3;
    float disc_inner = disc_radius - disc_width * 0.5;
    float disc_outer = disc_radius + disc_width * 0.5;
    vec3 disc_normal = vec3(0.0, 1.0, 0.0);
    float disc_thickness = 0.1;

    float center_dist = length(pos);
    float disc_dist = dot(disc_normal, pos);

    float rad_grad = 1.0 - clamp((center_dist - disc_inner) / disc_width * 0.5, 0.0, 1.0);
    float coverage = pcurve(rad_grad, 4.0, 0.9);
    disc_thickness *= rad_grad;

    vec3 dust_lit = vec3(1.0);
    vec3 dust_dark = vec3(0.0);

    float dust_glow = 1.0 / (pow(1.0 - rad_grad, 2.0) * 290.0 + 0.002);
    vec3 dust_color = dust_lit * dust_glow * 8.2;
    coverage = clamp(coverage * 0.7, 0.0, 1.0);

    float fade = pow((abs(center_dist - disc_inner) + 0.4), 4.0) * 0.04;
    float bloom_factor = 1.0 / (pow(disc_dist, 2.0) * 40.0 + fade + 0.000002);
    vec3 b = dust_lit * pow(bloom_factor, 1.5);

    b *= mix(vec3(1.7, 1.1, 1.0), vec3(0.5, 0.6, 1.0), vec3(pow(rad_grad, 2.0)));
    b *= mix(vec3(1.7, 0.5, 0.1), vec3(1.0, 1.0, 1.0), vec3(pow(rad_grad, 0.5)));

    dust_color = mix(dust_color, b * 150.0, clamp(1.0 - coverage * 1.0, 0.0, 1.0));
    coverage = clamp(coverage + bloom_factor * bloom_factor * 0.1, 0.0, 1.0);

    vec3 rad_coords = vec3(center_dist * 1.5 + 0.55,
                           atan2(-pos.x, -pos.z) * 1.5,
                           disc_dist * 1.5);
    rad_coords *= 0.95;

    float n1 = 1.0;
    vec3 rc = rad_coords + 0.0;             rc.y += time * SPEED;
    n1 *= noise(rc * 3.0) * 0.5 + 0.5;      rc.y -= time * SPEED;
    n1 *= noise(rc * 6.0) * 0.5 + 0.5;      rc.y += time * SPEED;
    n1 *= noise(rc * 12.0) * 0.5 + 0.5;     rc.y -= time * SPEED;
    n1 *= noise(rc * 24.0) * 0.5 + 0.5;     rc.y += time * SPEED;

    float n2 = 2.0;
    rc = rad_coords + 30.0;
    n2 *= noise(rc * 3.0) * 0.5 + 0.5;      rc.y += time * SPEED;
    n2 *= noise(rc * 6.0) * 0.5 + 0.5;      rc.y -= time * SPEED;
    n2 *= noise(rc * 12.0) * 0.5 + 0.5;     rc.y += time * SPEED;
    n2 *= noise(rc * 24.0) * 0.5 + 0.5;     rc.y -= time * SPEED;
    n2 *= noise(rc * 48.0) * 0.5 + 0.5;     rc.y += time * SPEED;
    n2 *= noise(rc * 92.0) * 0.5 + 0.5;     rc.y -= time * SPEED;

    dust_color *= n1 * 0.998 + 0.002;
    coverage *= n2;
    rad_coords.y += time * SPEED * 0.5;

    dust_color = max(vec3(0.0), dust_color);
    coverage *= 1200.0;
    coverage *= pcurve(rad_grad, 4.0, 0.9);

    return vec4(dust_color, coverage);

    // vec3 position = zeropos;
    // float lengthpos = length(position.xz);
    // float dist = min(0.25, lengthpos * (1.0 / SIZE) * 0.5) * SIZE * 0.4 * (1.0 / STEPS) / abs(ray.y);

    // position += dist * STEPS * ray * 0.5;


    // vec2 deltapos;
    // deltapos.x = -zeropos.z * 0.01 + zeropos.x;
    // deltapos.y = zeropos.x * 0.01 + zeropos.z;
    // deltapos = normalize(deltapos - zeropos.xz);

    // float parallel = dot(ray.xz, deltapos);
    // parallel /= sqrt(lengthpos);
    // parallel *= 0.7;
    // float redshift = parallel;// + 0.3;
    // redshift *= redshift;
    // redshift = clamp(redshift, 0.0, 1.0);

    // float dismix = clamp((lengthpos - SIZE * 2.0) * (1.0 / SIZE) * 0.24, 0.0, 1.0);
    // vec3 insidecol = mix(vec3(1.0, 0.8, 0.2), vec3(0.5, 0.13, 0.02) * 0.2, dismix);

    // insidecol *= mix(vec3(0.2, 0.2, 0.1), vec3(1.6, 2.4, 4.0), redshift);
    // insidecol *= 1.25;
    // redshift += 0.12;
    // redshift *= redshift;

    // vec4 o = vec4(0);
    // for (int i = 0; i < STEPS; i++) {
    //     position -= dist * ray;
    //     float intensity = clamp(1.0 - abs((i - 0.8) * (1.0 / STEPS) * 2.0), 0.0, 1.0);
    //     float lenpos = length(position.xz);
    //     float distmult = 1.0;
        
    //     distmult *= clamp((lenpos - SIZE * 0.75) * (1.0 / SIZE) * 1.5, 0.0, 1.0);
    //     distmult *= clamp((SIZE * 6.0 - lenpos) * (1.0 / SIZE) * 0.2, 0.0, 1.0);
    //     distmult *= distmult;

    //     float u = lenpos + time * SIZE * 0.3 + intensity * SIZE * 0.2;
    //     vec2 xy;
    //     float rot = mod(time * SPEED, 8192.0);
    //     float sinr = sin(rot);
    //     float cosr = cos(rot);
    //     xy.x = -position.z * sinr + position.x * cosr;
    //     xy.y = position.x * sinr + position.z * cosr;

    //     float x = abs(xy.x / xy.y);
    //     float angle = 0.02 * atan(x);

    //     vec3 rad_coords = vec3(lenpos * 1.5 + 0.55
    //                            atan2(-position.x, position.z) * 1.5
    //                            )
    //     float n1 = valnoise(vec2(angle, u * (1.0 / SIZE) * 0.05), f);
    //     n1 = n1 * 0.66 + 0.33 * valnoise(vec2(angle, u * (1.0 / SIZE) * 0.05), f * 2.0);

    //     float extrawidth = n1 * 1.0 * (1.0 - clamp(i * (1.0 / STEPS) * 2.0 - 1.0, 0.0, 1.0));
    //     float alpha = clamp(n1 * (intensity + extrawidth) * ((1.0 / SIZE) * 10.0 + 0.01) * dist * distmult, 0.0, 1.0);

    //     vec3 col = 2.0 * mix(vec3(0.3, 0.2, 0.15) * insidecol, insidecol, min(1.0, intensity * 2.0));
    //     o = clamp(vec4(col * alpha + o.rgb * (1.0 - alpha), o.a * (1.0 - alpha) + alpha), vec4(0.0), vec4(1.0));
    //     lenpos *= (1.0 / SIZE);
    //     o.rgb += redshift * (intensity * 1.0 + 0.5) * (1.0 / STEPS) * 100.0 * distmult / (lenpos * lenpos);
    // }
    // o.rgb = clamp(o.rgb - 0.005, 0.0, 1.0);
    // return o;
}

vec4 render(vec3 rd) {
    // vec2 angle = vec2(time * 0.1, 0.2);
    // rd = transform(rd, vec3(angle, 0.0), vec3(0));
    vec3 pos = camera_origin;

    vec4 col = vec4(0.0);
    vec4 glow = vec4(0.0);
    vec4 out_color = vec4(100.0);

    for (int disks = 0; disks < 20; disks++) {
        for (int h = 0; h < 6; h++) {
            float dotpos = dot(pos, pos);
            float invdist = inversesqrt(dotpos);
            float centdist = dotpos * invdist;
            float stepdist = 0.92 * abs(pos.y / rd.y);
            float farlim = centdist * 0.5;
            float closelim = centdist * 0.1 + 0.05 * centdist * centdist * (1.0 / SIZE);
            stepdist = min(stepdist, min(farlim, closelim));

            float invdist2 = invdist * invdist;
            float bendforce = stepdist * invdist2 * SIZE * 0.625;
            rd = normalize(rd - (bendforce * invdist) * pos);
            pos += stepdist * rd;
            glow += vec4(1.2, 1.1, 1.0, 1.0) * (0.01 * stepdist * invdist2 * invdist2 * clamp(centdist * 2.0 - 1.2, 0.0, 1.0));
        }

        float dist2 = length(pos);
        if (dist2 < SIZE * 0.1) {
            out_color = vec4(col.rgb * col.a + glow.rgb * (1.0 - col.a), 1.0);
            break;
        }
        else if (dist2 > SIZE * 1000.0) {
            vec4 bg = background(rd);
            out_color = vec4(col.rgb * col.a + bg.rgb * (1.0 - col.a) + glow.rgb * (1.0 - col.a), 1.0);
            break;
        }
        else if (abs(pos.y) <= SIZE * 0.002) {
            vec4 diskcol = raymarch_disk(rd, pos);
            pos.y = 0.0;
            pos += abs(SIZE * 0.001 / rd.y) * rd;
            col = vec4(diskcol.rgb * (1.0 - col.a) + col.rgb, col.a + diskcol.a * (1.0 - col.a));
        }
    }

    if (out_color.r == 100.0) out_color = vec4(col.rgb + glow.rgb * (col.a + glow.a), 1.0);
    
    col = out_color;
    // color grade
    col.rgb = pow(col.rgb, vec3(1.5));
    col.rgb = col.rgb / (1.0 + col.rgb);
    col.rgb = pow(col.rgb, vec3(1.0 / 1.5));

    col.rgb = mix(col.rgb, col.rgb * col.rgb * (3.0 - 2.0 * col.rgb), vec3(1.0));
    col.rgb = pow(col.rgb, vec3(1.3, 1.2, 1.0));
    col.rgb = clamp(col.rgb * 1.01, 0.0, 1.0);
    col.rgb = pow(col.rgb, vec3(0.7 / 2.2));
    return col;
}

void main() {
    vec2 uv = ((fragCoord * resolution.xy) - 0.5 * resolution.xy) / resolution.y;
    vec3 ray_direction = camera * normalize(vec3(uv, camera_focal_length));
    vec4 color = render(ray_direction);
    fragColor = color;
}