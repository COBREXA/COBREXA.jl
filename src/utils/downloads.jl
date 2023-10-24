
import SHA
import Downloads

function cached_download(
    url::AbstractString,
    path::AbstractString;
    output = missing,
    sha1 = missing,
    sha256 = missing,
    sha512 = missing,
    size = missing,
    kwargs...,
)
    ismissing(output) ||
        error("cached_download does not support the `output` keyword argument")

    ret = path

    if !isfile(path)
        ret = Downloads.download(url; output = path, kwargs...)
    end

    hash_matched = false

    function check_hash(hi, f, h, verifies = true)
        cksum = f(path)
        if cksum != h
            @warn "Downloaded file '$path' seems to be different from the expected one: got $hi '$cksum', expected $hi '$h'. This may be a cache corruption issue -- run `rm(\"$path\")` to flush the cache."
        end
        hash_checked |= verifies
    end

    compose(fs...) = foldl(ComposedFunction, fs)

    ismissing(sha1) || check_hash(:SHA1, compose(bytes2hex, SHA.sha1, open), sha1)
    ismissing(sha256) || check_hash(:SHA256, compose(bytes2hex, SHA.sha256, open), sha256)
    ismissing(sha512) || check_hash(:SHA512, compose(bytes2hex, SHA.sha512, open), sha512)
    ismissing(size) || check_hash(:size, filesize, size, false)

    hash_matched ||
        @warn "Downloaded file '$path' was not checked against any checksum. Add some to improve reproducibility."

    return ret
end
