@testitem "Upsampling" begin
    using DataInterpolations
    using Unitful

    # * Vector
    ts = (rand(100) .- 0.5) * 2π
    x = Timeseries(sinc.(ts), ts)
    @test_throws "`interpolate` only supports forward" TimeseriesTools.interpolate(x)

    sort!(ts)
    x = Timeseries(sinc.(ts), ts)
    itp = TimeseriesTools.interpolate(x)
    y = itp(dims(x) |> only)
    @test x ≈ y

    ts = (-π):0.1:π
    x = Timeseries(sinc.(ts), ts)
    itp = TimeseriesTools.interpolate(x)
    z = @test_nowarn upsample(x, 2)
    @test all(z[𝑡(At(ts))] .≈ x)
    @test z ≈ sinc.(lookup(z) |> only) atol = 1.0e-2

    # upsample is a thin wrapper over resample onto a `factor`× denser grid
    @test upsample(x, 2) == resample(x, upsample(dims(x, 𝑡), 2))
    @test step(lookup(upsample(x, 3), 𝑡)) ≈ step(ts) / 3

    # * Matrix
    x = Timeseries(randn(100, 100), 0.1:0.1:10, Var(1:100))
    y = @test_nowarn upsample(x, 2)
    z = @test_nowarn upsample(x, 2, dims = 1)
    @test y == z
    z = @test_nowarn upsample(x, 2, dims = (1, 2))
    @test all(y .≈ z[𝑡(At(lookup(y, 𝑡))), Var(At(lookup(y, Var)))])
    @test all(x .≈ z[𝑡(At(lookup(x, 𝑡))), Var(At(lookup(x, Var)))])
    @test DimensionalData.name.(dims(x)) == DimensionalData.name.(dims(z))
    @test length(dims(z, 1)) == length(dims(z, 2)) == 199

    # * 3D array
    x = Timeseries(randn(100, 100, 100), 0.1:0.1:10, Var(1:100), X(1:100))
    y = @test_nowarn upsample(x, 2)
    z = @test_nowarn upsample(x, 2, dims = 1)
    @test y == z
    z = @test_nowarn upsample(x, 2, dims = (3, 2, 1))
    @test all(y .≈ z[𝑡(At(lookup(y, 𝑡))), Var(At(lookup(y, Var))), X(At(lookup(y, X)))])
    @test all(x .≈ z[𝑡(At(lookup(x, 𝑡))), Var(At(lookup(x, Var))), X(At(lookup(x, X)))])
    @test DimensionalData.name.(dims(x)) == DimensionalData.name.(dims(z))
    @test length(dims(z, 1)) == length(dims(z, 2)) == 199

    # * Unitful data
    ts = ((-π):0.1:π) * u"s"
    x = Timeseries(sinc.(ustrip(ts)), ts) * u"V"
    itp = TimeseriesTools.interpolate(x)
    @test unit(eltype(itp(dims(x) |> only))) == NoUnits
    z = @test_nowarn upsample(x, 2)
    @test unit(eltype(z)) == u"V"
    @test unit(eltype(lookup(z) |> only)) == u"s"
    @test all(z[𝑡(At(ts))] .≈ x)
end

@testitem "Resampling" begin
    using DataInterpolations
    using Unitful

    # * Univariate, regular -> regular via Number dt
    ts = 0.0:0.1:10.0
    x = Timeseries(sin.(ts), ts)
    y = @test_nowarn resample(x, 0.25)
    @test y isa RegularTimeseries
    @test step(lookup(y, 𝑡)) ≈ 0.25
    @test first(lookup(y, 𝑡)) == first(ts)
    @test last(lookup(y, 𝑡)) ≤ last(ts)
    @test y ≈ sin.(lookup(y, 𝑡)) atol = 1.0e-3

    # * Number dt: grid stops at-or-before last sample
    y = resample(x, 0.3)
    @test last(lookup(y, 𝑡)) ≤ last(ts)
    @test last(lookup(y, 𝑡)) + 0.3 > last(ts)

    # * Univariate, irregular -> regular via Number dt
    using Random
    its = sort(rand(MersenneTwister(1), 200)) .* 10
    xi = Timeseries(sin.(its), its)
    yi = @test_nowarn resample(xi, 0.05)
    @test step(lookup(yi, 𝑡)) ≈ 0.05
    @test yi ≈ sin.(lookup(yi, 𝑡)) atol = 1.0e-2

    # * Univariate -> arbitrary target vector
    target = [0.1, 0.5, 1.0, 3.3, 7.7]
    yv = @test_nowarn resample(x, target)
    @test collect(lookup(yv, 𝑡)) == target
    @test yv ≈ sin.(target) atol = 1.0e-3

    # * Univariate -> target Dimension passed directly
    td = 𝑡(0.0:0.2:5.0)
    yd = @test_nowarn resample(x, td)
    @test lookup(yd, 𝑡) == lookup(td)
    @test yd ≈ sin.(lookup(yd, 𝑡)) atol = 1.0e-3

    # * Explicit interpolation type + forwarded kwarg
    yl = @test_nowarn resample(x, 0.25, LinearInterpolation)
    @test step(lookup(yl, 𝑡)) ≈ 0.25
    @test yl ≈ sin.(lookup(yl, 𝑡)) atol = 1.0e-2
    yc = @test_nowarn resample(
        x, 0.25, CubicSpline;
        extrapolation = ExtrapolationType.Constant
    )
    @test yc ≈ sin.(lookup(yc, 𝑡)) atol = 1.0e-3

    # * Multivariate (matrix)
    M = Timeseries(
        hcat(sin.(ts), cos.(ts), sinc.(ts)),
        ts,
        Var(1:3)
    )
    yM = @test_nowarn resample(M, 0.25)
    @test size(yM, 2) == 3
    @test step(lookup(yM, 𝑡)) ≈ 0.25
    @test lookup(yM, Var) == 1:3
    @test DimensionalData.name.(dims(yM)) == DimensionalData.name.(dims(M))
    @test yM[:, 1] ≈ sin.(lookup(yM, 𝑡)) atol = 1.0e-3
    @test yM[:, 2] ≈ cos.(lookup(yM, 𝑡)) atol = 1.0e-3

    # * 3D array
    T3 = Timeseries(
        randn(length(ts), 4, 5),
        ts, Var(1:4), X(1:5)
    )
    y3 = @test_nowarn resample(T3, 0.25)
    @test size(y3, 1) == length(lookup(y3, 𝑡))
    @test size(y3, 2) == 4 && size(y3, 3) == 5
    @test step(lookup(y3, 𝑡)) ≈ 0.25
    @test DimensionalData.name.(dims(y3)) == DimensionalData.name.(dims(T3))

    # * Unitful time dim (target Number with units)
    tsu = (0.0:0.1:10.0) * u"s"
    xu = Timeseries(sin.(ustrip.(tsu)), tsu)
    yu = @test_nowarn resample(xu, 0.25u"s")
    @test unit(eltype(lookup(yu, 𝑡))) == u"s"
    @test step(lookup(yu, 𝑡)) == 0.25u"s"
    @test ustripall(yu) ≈ sin.(ustrip.(lookup(yu, 𝑡))) atol = 1.0e-3

    # * Unitful values (and unitful time)
    xv = Timeseries(sin.(ustrip.(tsu)), tsu) * u"V"
    yv2 = @test_nowarn resample(xv, 0.25u"s")
    @test unit(eltype(lookup(yv2, 𝑡))) == u"s"
    @test unit(eltype(yv2)) == u"V"
    @test step(lookup(yv2, 𝑡)) == 0.25u"s"
    @test ustripall(yv2) ≈ sin.(ustrip.(lookup(yv2, 𝑡))) atol = 1.0e-3

    # * Unitful + target Dimension
    td_u = 𝑡((0.0:0.5:8.0) * u"s")
    yd_u = @test_nowarn resample(xu, td_u)
    @test lookup(yd_u, 𝑡) == lookup(td_u)
    @test ustripall(yd_u) ≈ sin.(ustrip.(lookup(yd_u, 𝑡))) atol = 1.0e-3
end

@testitem "Imputation" begin
    using DataInterpolations
    using Unitful

    ts = 0.0:0.1:10.0
    truth = sin.(ts)

    # * NaN values are filled
    v = collect(truth)
    bad = [10, 25, 60]
    v[bad] .= NaN
    x = Timeseries(v, ts)
    y = @test_nowarn impute(x)
    @test !any(isnan, parent(y))
    @test lookup(y, 𝑡) == lookup(x, 𝑡)
    @test y ≈ truth atol = 1.0e-2

    # * missing values are filled (Union eltype)
    vm = Vector{Union{Missing, Float64}}(truth)
    vm[bad] .= missing
    xm = Timeseries(vm, ts)
    ym = @test_nowarn impute(xm)
    @test !any(ismissing, parent(ym))
    @test eltype(ym) <: Real
    @test ym ≈ truth atol = 1.0e-2

    # * nothing values (Nothing in default replace set)
    vn = Vector{Any}(collect(truth))
    vn[bad] .= nothing
    xn = Timeseries(vn, ts)
    yn = @test_nowarn impute(xn)
    @test !any(isnothing, parent(yn))
    @test yn ≈ truth atol = 1.0e-2

    # * clean series passes through unchanged (to interpolation precision)
    xc = Timeseries(collect(truth), ts)
    yc = @test_nowarn impute(xc)
    @test yc ≈ xc atol = 1.0e-6
    @test lookup(yc, 𝑡) == lookup(xc, 𝑡)

    # * custom replace: a sentinel value
    vs = collect(truth)
    vs[bad] .= -999.0
    xs = Timeseries(vs, ts)
    ys = @test_nowarn impute(xs; replace = [-999.0])
    @test all(!=(-999.0), parent(ys))
    @test ys ≈ truth atol = 1.0e-2

    # * custom replace: a type (Complex)
    vcpx = Vector{Any}(collect(truth))
    vcpx[bad] = [1.0 + 2.0im for _ in bad]
    xcpx = Timeseries(vcpx, ts)
    ycpx = @test_nowarn impute(xcpx; replace = [NaN, Nothing, Missing, Complex])
    @test !any(z -> z isa Complex, parent(ycpx))
    @test ycpx ≈ truth atol = 1.0e-2

    # * explicit interpolation type forwarded
    yl = @test_nowarn impute(x, LinearInterpolation)
    @test !any(isnan, parent(yl))
    @test yl ≈ truth atol = 5.0e-2

    # * multivariate: each column imputed along dims = 1
    M = hcat(collect(truth), collect(cos.(ts)))
    M[bad, 1] .= NaN
    M[[15, 40], 2] .= NaN
    xM = Timeseries(M, ts, Var(1:2))
    yM = @test_nowarn impute(xM)
    @test !any(isnan, parent(yM))
    @test size(yM) == size(xM)
    @test DimensionalData.name.(dims(yM)) == DimensionalData.name.(dims(xM))
    @test yM[:, 1] ≈ truth atol = 1.0e-2
    @test yM[:, 2] ≈ cos.(ts) atol = 1.0e-2

    # * unitful values and unitful time
    tsu = (0.0:0.1:10.0) * u"s"
    vu = collect(sin.(ustrip.(tsu)))
    vu[bad] .= NaN
    xu = Timeseries(vu, tsu) * u"V"
    yu = @test_nowarn impute(xu)
    @test unit(eltype(yu)) == u"V"
    @test unit(eltype(lookup(yu, 𝑡))) == u"s"
    @test !any(z -> isnan(ustrip(z)), parent(yu))
    @test ustripall(yu) ≈ sin.(ustrip.(tsu)) atol = 1.0e-2
end

@testitem "Downsampling" begin
    using DSP
    using FFTW
    using Statistics
    using Unitful

    fs = 200.0
    t = 0.0:(1 / fs):2.0   # an AbstractRange, so a RegularTimeseries

    # * Low-frequency signal: well below the new Nyquist, both modes preserve it
    lo = Timeseries(sin.(2π * 5 .* t), t)
    @test lo isa RegularTimeseries
    d = @test_nowarn downsample(lo, 4)
    @test d isa RegularTimeseries
    @test samplingrate(d) ≈ fs / 4
    @test step(lookup(d, 𝑡)) ≈ 4 / fs
    @test first(lookup(d, 𝑡)) == first(t)
    # compare on the interior, away from FIR filter edge transients
    interior = 6:(length(d) - 5)
    @test d[interior] ≈ sin.(2π * 5 .* lookup(d, 𝑡))[interior] atol = 1.0e-2

    # * factor = 1 is a no-op
    @test downsample(lo, 1) === lo
    @test_throws ArgumentError downsample(lo, 0)

    # * antialias = false is exact plain decimation
    d0 = @test_nowarn downsample(lo, 4; antialias = false)
    @test parent(d0) == parent(lo)[1:4:end]
    @test step(lookup(d0, 𝑡)) ≈ 4 / fs

    # * Aliasing: a 70 Hz tone, decimating by 4 (new Nyquist 25 Hz) must fold without filter
    hi = Timeseries(sin.(2π * 70 .* t), t)
    naive = downsample(hi, 4; antialias = false)
    filt = downsample(hi, 4; antialias = true)
    # the filtered version suppresses the out-of-band tone -> much less power.
    # (loose factor leaves margin for FIR transition-band leakage and edge transients)
    @test var(parent(filt)) < 0.5 * var(parent(naive))
    # naive decimation retains near-full power: the 70 Hz tone aliased fully into the
    # kept band, so its variance stays close to a unit sine's ~0.5.
    @test var(parent(naive)) > 0.3

    # * Multivariate: each column downsampled along the time dim
    M = Timeseries(hcat(sin.(2π * 5 .* t), sin.(2π * 3 .* t)), t, Var(1:2))
    dM = @test_nowarn downsample(M, 4)
    @test size(dM, 2) == 2
    @test samplingrate(dM) ≈ fs / 4
    @test DimensionalData.name.(dims(dM)) == DimensionalData.name.(dims(M))
    intM = 6:(size(dM, 1) - 5)
    @test dM[intM, 1] ≈ sin.(2π * 5 .* lookup(dM, 𝑡))[intM] atol = 1.0e-2
    @test dM[intM, 2] ≈ sin.(2π * 3 .* lookup(dM, 𝑡))[intM] atol = 1.0e-2

    # * Unitful values and time
    xu = Timeseries(sin.(2π * 5 .* t), t * u"s") * u"V"
    du = @test_nowarn downsample(xu, 4)
    @test unit(eltype(du)) == u"V"
    @test unit(eltype(lookup(du, 𝑡))) == u"s"
    @test step(lookup(du, 𝑡)) == (4 / fs) * u"s"
end

@testitem "DataInterpolationsNDExt" begin
    using DataInterpolationsND
    using Unitful

    tx = 0.0:0.5:5.0
    ty = 0.0:0.5:4.0
    # * Bilinear field: linear ND interpolation is exact on it
    bilin(x, y) = 1.0 + 2.0x - 0.5y + 0.3x * y
    U = [bilin(x, y) for x in tx, y in ty]
    X = Timeseries(U, 𝑡(tx), Var(ty))

    # * interpolate -> NDInterpolation, callable on target dims
    itp = @test_nowarn TimeseriesTools.interpolate(X, LinearInterpolationDimension)
    @test itp isa NDInterpolation
    # evaluate on the original grid -> recovers the data
    Y = itp((𝑡(tx), Var(ty)))
    @test Y ≈ U atol = 1.0e-10
    # evaluate on a finer grid -> matches the analytic field (bilinear is exact)
    fxt = 𝑡(0.0:0.25:5.0)
    fyt = Var(0.0:0.25:4.0)
    Yf = itp((fxt, fyt))
    @test Yf ≈ [bilin(x, y) for x in lookup(fxt), y in lookup(fyt)] atol = 1.0e-10
    @test DimensionalData.name.(dims(Yf)) == DimensionalData.name.(dims(X))

    # * upsample, joint over both axes
    Z = @test_nowarn upsample(X, 2, LinearInterpolationDimension)
    @test length(dims(Z, 1)) == 2 * length(tx) - 1
    @test length(dims(Z, 2)) == 2 * length(ty) - 1
    @test DimensionalData.name.(dims(Z)) == DimensionalData.name.(dims(X))
    @test Z ≈ [bilin(x, y) for x in lookup(Z, 1), y in lookup(Z, 2)] atol = 1.0e-10
    # original samples preserved
    @test Z[𝑡(At(collect(tx))), Var(At(collect(ty)))] ≈ U atol = 1.0e-10

    # * upsample along a single axis only
    Z1 = @test_nowarn upsample(X, 2, LinearInterpolationDimension; dims = 1)
    @test length(dims(Z1, 1)) == 2 * length(tx) - 1
    @test length(dims(Z1, 2)) == length(ty)

    # * per-axis dimension specs (tuple)
    itp2 = @test_nowarn TimeseriesTools.interpolate(
        X, (LinearInterpolationDimension, LinearInterpolationDimension)
    )
    @test itp2 isa NDInterpolation

    # * resample onto explicit target lookups
    R = @test_nowarn resample(X, (0.0:1.0:5.0, 0.0:1.0:4.0), LinearInterpolationDimension)
    @test lookup(R, 1) == 0.0:1.0:5.0
    @test lookup(R, 2) == 0.0:1.0:4.0
    @test R ≈ [bilin(x, y) for x in lookup(R, 1), y in lookup(R, 2)] atol = 1.0e-10

    # * resample with a single Number (common period on all axes)
    Rn = @test_nowarn resample(X, 0.25, LinearInterpolationDimension)
    @test step(lookup(Rn, 1)) ≈ 0.25
    @test step(lookup(Rn, 2)) ≈ 0.25
    @test last(lookup(Rn, 1)) ≤ last(tx)
    @test last(lookup(Rn, 2)) ≤ last(ty)
    @test Rn ≈ [bilin(x, y) for x in lookup(Rn, 1), y in lookup(Rn, 2)] atol = 1.0e-10

    # * unitful values and unitful lookups
    Xu = Timeseries(U, 𝑡(tx * u"s"), Var(ty)) * u"V"
    Zu = @test_nowarn upsample(Xu, 2, LinearInterpolationDimension)
    @test unit(eltype(Zu)) == u"V"
    @test unit(eltype(lookup(Zu, 1))) == u"s"
    @test ustripall(Zu) ≈ [bilin(x, y) for x in ustrip.(lookup(Zu, 1)), y in lookup(Zu, 2)] atol = 1.0e-8

    # * ConstantInterpolationDimension also works as a grid interpolator
    Rc = @test_nowarn resample(X, (0.0:1.0:5.0, 0.0:1.0:4.0), ConstantInterpolationDimension)
    @test lookup(Rc, 1) == 0.0:1.0:5.0
    @test lookup(Rc, 2) == 0.0:1.0:4.0
    # constant interpolation reproduces grid points exactly
    @test Rc[𝑡(At(0.0)), Var(At(0.0))] ≈ bilin(0.0, 0.0) atol = 1.0e-10

    # * BSplineInterpolationDimension: degree 1 is exact piecewise-linear interpolation
    itpb1 = @test_nowarn TimeseriesTools.interpolate(X, BSplineInterpolationDimension; degree = 1)
    @test itpb1 isa NDInterpolation
    Yb1 = itpb1((𝑡(tx), Var(ty)))
    @test Yb1 ≈ U atol = 1.0e-10                      # degree 1 -> passes through samples
    Ub1 = @test_nowarn upsample(X, 2, BSplineInterpolationDimension; degree = 1)
    @test Ub1 ≈ [bilin(x, y) for x in lookup(Ub1, 1), y in lookup(Ub1, 2)] atol = 1.0e-10

    # * BSplineInterpolationDimension: degree 2 runs (smoothing, not exact) and stays bounded
    Yb2 = @test_nowarn upsample(X, 2, BSplineInterpolationDimension; degree = 2)
    @test length(dims(Yb2, 1)) == 2 * length(tx) - 1
    @test DimensionalData.name.(dims(Yb2)) == DimensionalData.name.(dims(X))
    @test all(minimum(U) - 1 .≤ Yb2 .≤ maximum(U) + 1)   # bounded, smoothed
end
