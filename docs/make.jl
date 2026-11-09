using Documenter, DocumenterCitations, DocumenterInterLinks
using StandardHaloModelsPP

# Documenter Meta & Plugins ========================================================
DocMeta.setdocmeta!(StandardHaloModelsPP, :DocTestSetup, :(using StandardHaloModelsPP); recursive=true)

makebib() = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

link = InterLinks(
    "Unitful" => (
        "https://juliaphysics.github.io/Unitful.jl/stable",
        "https://juliaphysics.github.io/Unitful.jl/stable/objects.inv"
    ),
#    "NaturallyUnitful" => (
#        "https://juliaphysics.github.io/NaturallyUnitful.jl/stable",
#        "https://juliaphysics.github.io/NaturallyUnitful.jl/stable/objects.inv"
#    ),
    "Measurements" => (
        "https://juliaphysics.github.io/Measurements.jl/stable/",
        "https://juliaphysics.github.io/Measurements.jl/stable/objects.inv"
    ),
    "Integrals" => (
        "https://docs.sciml.ai/Integrals/stable/",
        "https://docs.sciml.ai/Integrals/stable/objects.inv"
    )
)


# Asymptote Figures ================================================================
const FIG_DIR = joinpath(@__DIR__, "src", "figs")
const FIG_SRC = ["gal.asy", "eph.asy"]
const FIG_FORMATS = [
    "html" => ["-render", "4", "-offline"],
    "pdf"  => ["-render", "0", "-tex", "pdflatex"],
]

const ASY_TOOLCHAIN_HINT = """
Building the StandardHaloModelsPP.jl documentation requires Asymptote and a LaTeX
installation on PATH.

Set the ASY environment variable to override `asy` executable discovery, and
SHM_REBUILD_FIGS=true to force regeneration of up-to-date figures.
"""

figsrc(src) = joinpath(FIG_DIR, src)
figout(src, fmt) = joinpath(FIG_DIR, first(splitext(src)) * "." * fmt)

function asy_toolchain(formats)
    asy = Sys.which(get(ENV, "ASY", "asy"))
    needed = ["latex", "dvips", "gs"]
    "pdf" in formats && push!(needed, "pdflatex")
    absent = [t for t in needed if Sys.which(t) === nothing]
    asy === nothing && pushfirst!(absent, "asy")

    isempty(absent) || error("""
        Cannot build Asymptote doc figures - missing: $(join(absent, ", ")).

        $ASY_TOOLCHAIN_HINT""")

    ver = IOBuffer()
    run(pipeline(ignorestatus(`$asy -version`); stdout = ver, stderr = ver))
    @info "Asymptote toolchain" asy version = first(split(String(take!(ver)), '\n'))
    return asy
end

function asy_render(asy::AbstractString, src::AbstractString,
                    fmt::AbstractString, flags::Vector{String})
    cmd = Cmd(`$asy -noV -safe -f $fmt $flags $src`; dir = FIG_DIR)
    @info "Rendering Asymptote figure" src fmt
    log = IOBuffer()
    try
        run(pipeline(cmd; stdout = log, stderr = log))
    catch err
        error("""
            Asymptote failed rendering $src as $fmt ($err).
              command: $cmd
              output:
            $(String(take!(log)))""")
    end
    out = String(take!(log))
    isfile(figout(src, fmt)) ||
        error("asy exited cleanly but produced no $(figout(src, fmt)); log:\n$out")
    isempty(strip(out)) ||
        @warn "Asymptote reported diagnostics - check the figure's labels" src fmt log = out
    return nothing
end

function asy_build(; formats = FIG_FORMATS,
                     force::Bool = get(ENV, "SHM_REBUILD_FIGS", "false") == "true")
    jobs = [(src, fmt, flags) for src in FIG_SRC for (fmt, flags) in formats
            if force || !isfile(figout(src, fmt)) ||
               mtime(figsrc(src)) > mtime(figout(src, fmt))]

    if isempty(jobs)
        @info "Asymptote figures up to date"
        return nothing
    end

    asy = asy_toolchain(first.(formats))
    foreach(job -> asy_render(asy, job...), jobs)
    @info "Rendered $(length(jobs)) Asymptote figure(s)"
    return nothing
end

asy_build()

# Make Documentation! ==============================================================

const PAGES = [
    "index.md",
    "ephemeris.md",
    "integrals.md",
    "performance.md",
    "Reference" => "ref.md",
]

# PDF
makedocs(;
    modules=[StandardHaloModelsPP],
    authors="Thomas (Tom) C. Gorordo",
    sitename="StandardHaloModelsPP.jl",
    build="pdf",
    format=Documenter.LaTeX(show_log = true),
    pages=PAGES,
    pagesonly = true,
    plugins = [makebib(), link]
)

# HTML/Site
makedocs(;
    modules=[StandardHaloModelsPP],
    authors="Thomas (Tom) C. Gorordo",
    sitename="StandardHaloModelsPP.jl",
    build="html",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
        prettyurls=true
    ),
    pages=PAGES,
    pagesonly = true,
    plugins = [makebib(), link]
)

let pdfs = filter(endswith(".pdf"), readdir(joinpath(@__DIR__, "pdf"); join = true))
    for p in pdfs
        cp(p, joinpath(@__DIR__, "html", basename(p)); force = true)
    end
end

deploydocs(repo="github.com/tgorordo/StandardHaloModelsPP.jl.git", target="html")



