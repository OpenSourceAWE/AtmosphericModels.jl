# SPDX-FileCopyrightText: 2026 Uwe Fechner
# SPDX-License-Identifier: MIT

using Test

"""
    commit_release(repo, project_version, changelog_version, notes)

Commit a Project.toml and CHANGELOG.md for one release to `repo` and return its sha.
"""
function commit_release(repo, project_version, changelog_version, notes)
    write(joinpath(repo, "Project.toml"), "version = \"$project_version\"\n")
    write(joinpath(repo, "CHANGELOG.md"),
          "## AtmosphericModels v$changelog_version 2026-11-08\n\n- $notes\n")
    run(`git -C $repo add -A`)
    identity = `-c user.name=test -c user.email=test@example.com`
    run(`git -C $repo $identity commit -q -m release`)
    return readchomp(`git -C $repo rev-parse HEAD`)
end

"""
    run_release_from_temporary_copy(project_version, changelog_version)

Run a temporary copy of `bin/release` for a commit that is neither checked out nor fetched
in the clone it runs in, with `gh` faked. Return its exit code, its output, the commit's sha
and the arguments `gh` got (`nothing` if it was not called).
"""
function run_release_from_temporary_copy(project_version, changelog_version)
    mktempdir() do box
        origin, clone, script_dir, bin_dir =
            joinpath.(box, ("origin", "clone", "script", "bin"))
        mkpath.((origin, script_dir, bin_dir))
        run(`git -C $origin init -q -b main`)
        sha = commit_release(origin, project_version, changelog_version,
                             "notes of this release")
        commit_release(origin, "9.9.9", "9.9.9", "notes of a later release")
        run(`git clone -q --depth 1 file://$origin $clone`)

        gh_log = joinpath(box, "gh.log")
        write(joinpath(bin_dir, "gh"), """
              #!/bin/bash
              printf '%s\\n' "\$@" > $gh_log
              echo https://example.com/posted
              """)
        chmod(joinpath(bin_dir, "gh"), 0o755)
        release = cp(joinpath(@__DIR__, "..", "bin", "release"),
                     joinpath(script_dir, "release"))

        cmd = addenv(`$release --yes --sha $sha`, "PATH" => "$bin_dir:$(ENV["PATH"])",
                     "GH_REPO" => "OpenSourceAWE/AtmosphericModels.jl")
        output = IOBuffer()
        process = run(pipeline(ignorestatus(Cmd(cmd; dir=clone));
                               stdout=output, stderr=output))
        gh_args = isfile(gh_log) ? read(gh_log, String) : nothing
        return (; exit_code=process.exitcode, output=String(take!(output)), sha, gh_args)
    end
end

@testset "bin/release registers the given commit from a temporary copy of the script" begin
    result = run_release_from_temporary_copy("0.3.13", "0.3.13")
    @test result.exit_code == 0
    @test contains(result.output, "Posted: https://example.com/posted")
    @test contains(result.gh_args, "repos/{owner}/{repo}/commits/$(result.sha)/comments")
    @test contains(result.gh_args, "body=@JuliaRegistrator register")
    @test contains(result.gh_args, "notes of this release")
    @test !contains(result.gh_args, "notes of a later release")
end

@testset "bin/release refuses a commit whose changelog and Project.toml disagree" begin
    result = run_release_from_temporary_copy("0.3.13", "0.3.12")
    @test result.exit_code != 0
    @test contains(result.output, "Version mismatch")
    @test result.gh_args === nothing
end
