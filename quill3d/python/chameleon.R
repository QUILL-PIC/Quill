dyn.load("libchameleon-R.so")

# see for writing of R bindings:
# http://adv-r.had.co.nz/C-interface.html

foo <- function(i) {
    print("chameleon: foo() loaded:")
    print(is.loaded("foo"))
    return(.Call("foo", i = as.integer(i)))
}

configure <- function(log, flag = FALSE) {
    .Call("configure_fr", log = as.character(log), flag = as.logical(flag))
}

get <- function(name) {
    .Call("get_fr", name = as.character(name))
}

read2d <- function(filename, plane) {
    xname <- substring(plane, 1, 1)
    yname <- substring(plane, 2, 2)
    dx <- get(paste("d", xname, sep = ""))
    dy <- get(paste("d", yname, sep = ""))
    nx <- as.integer(get(paste("n", xname, sep = "")))
    ny <- as.integer(get(paste("n", yname, sep = "")))
    xl <- dx * (nx - 1)
    yl <- dy * (ny - 1)
    x <- seq(0, xl, len = nx)
    y <- seq(0, yl, len = ny)
    a <- expand.grid(y = y, x = x)
    a$z <- .Call("read_fr", filename = as.character(filename), plane =
                 as.character(plane))
    return(a);
}
