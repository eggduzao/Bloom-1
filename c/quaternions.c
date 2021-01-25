#include <stdio.h>
#ifdef __arm__
#define LIBQUAT_ARM
#include <arm_math.h>
#else
#include <math.h>
#endif
#include <assert.h>

#include "geometry_datatypes.h"

quaternion quat_from_euler(const euler *in_euler)
{
    quaternion out_quat;

    if (!in_euler)
    {
        fprintf(stderr, "%s: Invalid input.", __func__);
        return out_quat;
    }

    quaternion temp;

    float cy = cosf(in_euler->yaw * 0.5f);
    float sy = sinf(in_euler->yaw * 0.5f);
    float cp = cosf(in_euler->pitch * 0.5f);
    float sp = sinf(in_euler->pitch * 0.5f);
    float cr = cosf(in_euler->roll * 0.5f);
    float sr = sinf(in_euler->roll * 0.5f);

    temp.w = cr * cp * cy + sr * sp * sy;
    temp.q1 = sr * cp * cy - cr * sp * sy;
    temp.q2 = cr * sp * cy + sr * cp * sy;
    temp.q3 = cr * cp * sy - sr * sp * cy;

    return temp;
}

euler euler_from_quat(const quaternion *in_quat)
{
    euler out_euler;
    if (!in_quat)
    {
        fprintf(stderr, "%s: Invalid input.", __func__);
        return out_euler;
    }

    out_euler.roll = atan2f(
        2.f * (in_quat->w * in_quat->q1 + in_quat->q2 * in_quat->q3),
        1.f - 2.f * (in_quat->q1 * in_quat->q1 + in_quat->q2 * in_quat->q2));
    out_euler.pitch =
        asinf(2.f * (in_quat->w * in_quat->q2 + in_quat->q1 * in_quat->q3));
    out_euler.yaw = atan2f(
        2.f * (in_quat->w * in_quat->q3 + in_quat->q1 * in_quat->q2),
        1.f - 2.f * (in_quat->q2 * in_quat->q2 + in_quat->q3 * in_quat->q3));

    return out_euler;
}

quaternion quaternion_multiply(const quaternion *in_quat1,
                               const quaternion *in_quat2)
{
    quaternion out_quat;
    if (!in_quat1 || !in_quat2)
    {
        fprintf(stderr, "%s: Invalid input.", __func__);
        return out_quat;
    }

    out_quat.w = in_quat1->w * in_quat2->w - in_quat1->q1 * in_quat2->q1 -
                 in_quat1->q2 * in_quat2->q2 - in_quat1->q3 * in_quat2->q3;
    out_quat.q1 = in_quat1->w * in_quat2->q1 + in_quat1->q1 * in_quat2->w +
                  in_quat1->q2 * in_quat2->q3 - in_quat1->q3 * in_quat2->q2;
    out_quat.q2 = in_quat1->w * in_quat2->q2 - in_quat1->q1 * in_quat2->q3 +
                  in_quat1->q2 * in_quat2->w + in_quat1->q3 * in_quat2->q1;
    out_quat.q3 = in_quat1->w * in_quat2->q3 + in_quat1->q1 * in_quat2->q2 -
                  in_quat1->q2 * in_quat2->q1 + in_quat1->q3 * in_quat2->w;

    return out_quat;
}

