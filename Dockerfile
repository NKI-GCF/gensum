FROM rust:latest as buildenv

WORKDIR /usr/app/src
COPY ./ /usr/app/src

RUN apt-get update && apt-get -y install cmake build-essential clang && \
    rm -rf /var/lib/apt/lists/* && \
    rustup component add rustfmt 

RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/rust/target \
    cargo build --release

FROM debian:bullseye-slim as runner
WORKDIR /root
COPY --from=buildenv /usr/app/src/target/release/ /usr/local/bin/
RUN chmod +x /usr/local/bin/gensum
CMD /usr/local/bin/gensum
