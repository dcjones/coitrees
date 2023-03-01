use std::arch::aarch64::*;

use std::mem::transmute;
fn count_bits(bits: uint32x4_t) -> usize {
    unsafe {
        let t2 = vreinterpretq_u8_u32(bits);
        let t3 = vcntq_u8(t2);
        let sum = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(t3)));
        vgetq_lane_u64(sum, 0) as usize + vgetq_lane_u64(sum, 1) as usize
    }
}

fn main() {
    let first_data = [1, 9, 3, 4];
    let second_data = [5, 6, 7, 8];

    unsafe {
        let firsts = vld1q_s32(first_data.as_ptr());
        let seconds = vld1q_s32(second_data.as_ptr());

        let value: u128 = transmute(firsts);
        println!("value: {:#0128b}", value);

        print_4_i32(firsts);
        print_4_i32(seconds);

        let cmp = vcgtq_s32(firsts, seconds);

        println!("len: {}", count_bits(cmp));
        print_4_u32(cmp);

        let naw = vmovn_u32(cmp);

        print_4_u16(naw);

        let value: u64 = transmute(naw);

        let mask: u64 = 0xffff;
        println!("mask : {:#b}", value & mask);
        println!("mask  2: {:#064b}", value & mask << 16);
        println!("mask  3: {:#b}", value & mask << 32);
        println!("mask  4: {:#b}", value & mask << 48);

        // println!("mask : {:#b}", value & mask);

        // eprintln!("DEBUGPRINT[5]: neon.rs:25: value={:#b}", value);

        // let a = "11111111111111111111111111111111";

        // let cmp = vreinterpretq_u8_u32(cmp);
        // let cmp = vaddvq_u8(cmp);
    }
}

fn print_4_u32(cmp: uint32x4_t) {
    unsafe {
        let number1 = vgetq_lane_u32(cmp, 0);
        eprintln!("DEBUGPRINT[1]: neon.rs:12: number1={:#b}", number1);
        let number2 = vgetq_lane_u32(cmp, 1);
        eprintln!("DEBUGPRINT[2]: neon.rs:14: number2={:#b}", number2);
        let number3 = vgetq_lane_u32(cmp, 2);
        eprintln!("DEBUGPRINT[3]: neon.rs:16: number3={:#b}", number3);
        let number4 = vgetq_lane_u32(cmp, 3);
        eprintln!("DEBUGPRINT[4]: neon.rs:18: number4={:#b}", number4);
    }
}

fn print_4_i32(cmp: int32x4_t) {
    unsafe {
        let number1 = vgetq_lane_s32(cmp, 0);
        eprintln!("DEBUGPRINT[1]: neon.rs:12: number1={:#b}", number1);
        let number2 = vgetq_lane_s32(cmp, 1);
        eprintln!("DEBUGPRINT[2]: neon.rs:14: number2={:#b}", number2);
        let number3 = vgetq_lane_s32(cmp, 2);
        eprintln!("DEBUGPRINT[3]: neon.rs:16: number3={:#b}", number3);
        let number4 = vgetq_lane_s32(cmp, 3);
        eprintln!("DEBUGPRINT[4]: neon.rs:18: number4={:#b}", number4);
    }
}

fn print_4_u16(cmp: uint16x4_t) {
    unsafe {
        let number1 = vget_lane_u16(cmp, 0);
        eprintln!("DEBUGPRINT[1]: neon.rs:12: number1={:#b}", number1);
        let number2 = vget_lane_u16(cmp, 1);
        eprintln!("DEBUGPRINT[2]: neon.rs:14: number2={:#b}", number2);
        let number3 = vget_lane_u16(cmp, 2);
        eprintln!("DEBUGPRINT[3]: neon.rs:16: number3={:#b}", number3);
        let number4 = vget_lane_u16(cmp, 3);
        eprintln!("DEBUGPRINT[4]: neon.rs:18: number4={:#b}", number4);
    }
}
