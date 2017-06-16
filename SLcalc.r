## matrix & vector math for simple ray-trace problem
## for a bounded emitting plane, calculate the illuminance on a receiver plane

## Helper function from RSEIS (3D vector cross product)
xprod <- function(v1, v2) {
    x = v1[2] * v2[3] - v1[3] * v2[2]
    y = v1[3] * v2[1] - v1[1] * v2[3]
    z = v1[1] * v2[2] - v2[1] * v1[2]
    return(c(x, y, z))
}

### Helper function returns the length of a 3D vector
vlength <- function(v) {
    return(sqrt(sum(v*v)))
}

## Helper function returns the unit vector joining two points in 3D space
lvector <- function(r1, r2) {
    l = vlength(r2-r1)
    return((r2-r1)/l)
}

## Helper function returns unit normal of the plane defined by two vectors in the plane
pnormal <- function(v1, v2) {
    l = vlength(xprod(v1, v2))
    x = xprod(v1, v2)[1]/l
    y = xprod(v1, v2)[2]/l
    z = xprod(v1, v2)[3]/l
    return(c(x, y, z))
}

## Helper function returns cosine of the intersection angle between a line 
## and a plane (defined by their direction and normal vectors, respectively)
plcos <- function(n, u){
    n = n/vlength(n)
    u = u/vlength(u)
    return(sum(n*u))
}

## Helper function rotates a 3D vector about an axis
## will make it convenient to replicate rotationally symmetric objects

## Calculate vertices (as 3D vectors) of a skylight type entity
## based on opening width, sidewall depth, & sidewall angle in degrees
## assumed to be parallel to the x-y plane, centered on the z-axis,
## and the center of the opening is (0,0,0)
corners <- function (w, d, alpha = 0){
    M <- matrix(nrow=8, ncol=3)
    s = w - 2*d*sin(alpha*pi/180)
    M[1,] = c(w/2, w/2, 0)
    M[2,] = c(w/2, -w/2, 0)
    M[3,] = c(-w/2, -w/2, 0)
    M[4,] = c(-w/2, w/2, 0)
    M[5,] = c(s/2, s/2, d)
    M[6,] = c(s/2, -s/2, d)
    M[7,] = c(-s/2, -s/2, d)
    M[8,] = c(-s/2, s/2, d)
    colnames(M)<-c("x", "y", "z")
    rownames(M)<-c("B1", "B2", "B3", "B4","T1", "T2", "T3", "T4")
    return(M)
}

 
## Function that defines a bounded plane based on the 4 vertices

